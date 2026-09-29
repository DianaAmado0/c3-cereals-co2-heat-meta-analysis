library(readxl)
library(metafor)
library(writexl)

if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# import dataset
source_file <- "dataset.xlsx"

# folder to save files
base_dir    <- "RESULTS"

N_GRID  <- 100
DIGITS  <- 6

SPECIE_COL <- "genotype"
SPECIE_PARENT <- "crop_type"
MODS_CAT <- c("experiment_type", "crop_type", SPECIE_COL)
MIN_STUDIES_LEVEL <- 2
MIN_LEVELS        <- 2

MAKE_ASSOC        <- TRUE
ASSOC_MIN_UNITS   <- 5
ASSOC_MIN_STUDIES <- 4
ASSOC_P_ADJUST    <- "BH"
ASSOC_MAX_PAIRS   <- 40

# Each trait pair is tested only ONCE (not in both directions), with the
# predictor chosen by causal order below. Testing both directions would
# double-count the same relationship in the FDR correction.
ASSOC_BOTH_DIRECTIONS <- FALSE
ASSOC_PREDICTOR_ORDER <- c("^grain yield|^yield", "biomass", "harvest index",
                           "photosynth", "stomatal", "thousand grain|grain weight",
                           "starch")

MAKE_REPEATED_GENO <- TRUE
GENO_MIN_STUDIES   <- 2

MAKE_IONOME_BY_GENO <- TRUE
IONOME_CATEGORY     <- "Grain ionome"
IONOME_MIN_OBS      <- 3
IONOME_MIN_STUDIES  <- 2

YIELD_PATTERN       <- "yield|rendimiento|grain mass|grain weight|seed weight|biomass"
DILUTION_CATEGORIES <- c("Grain ionome", "Grain quality")

YEAR_CANDIDATES <- c("year", "Year", "publication_year", "pub_year",
                     "anio", "ano", "yr")

clean_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))

signif_stars <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01,  "**",
                       ifelse(p < 0.05,  "*",
                              ifelse(p < 0.1,   ".", "ns")))))
}

to_numeric <- function(x) {
  if (is.numeric(x)) return(x)
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "na", "n/a", "n/d", "-", "--", ".")] <- NA
  x <- gsub("\\s", "", x)
  x <- gsub(",", ".", x, fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

bind <- function(l) if (length(l) > 0) do.call(rbind, l) else data.frame()

pct <- function(x) (exp(x) - 1) * 100

# unit key: same experimental condition = comparable across traits
make_unit_key <- function(d) {
  paste(d$source_file, d$crop_type, d[[SPECIE_COL]],
        d$experiment_type, round(d$delta_co2_ppm), sep = "|")
}

# ---- heterogeneity (I2), partitioned into between/within studies ----
i2_multilevel <- function(m, vi) {
  out <- list(I2_total = NA, I2_between = NA, I2_within = NA)
  res <- tryCatch({
    W  <- diag(1 / vi)                # weights = 1/variance
    X  <- model.matrix(m)             # design matrix
    P  <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
    s2 <- (m$k - m$p) / sum(diag(P))  # typical sampling variance
    tot <- sum(m$sigma2)              # sum of variance components
    list(I2_total   = 100 * tot / (tot + s2),
         I2_between = 100 * m$sigma2[1] / (tot + s2),
         I2_within  = if (length(m$sigma2) > 1) 100 * m$sigma2[2] / (tot + s2) else NA)
  }, error = function(e) NULL)
  if (!is.null(res)) out <- res
  out
}

# Multilevel is used to avoid pseudoreplication and to account for the
# dependence among outcomes that come from the same study.
# random = ~ 1 | study_id/obs_id  adds random effects at both levels
# (between studies and within studies) (Credé et al. 2010).
# yi   = observed effect size estimates
# V    = corresponding sampling variances
# REML = restricted maximum likelihood (unbiased heterogeneity estimates)

fit_trait <- function(sub) {
  sub$study_id <- factor(sub$study_number)    # identification of each study
  sub$obs_id   <- factor(seq_len(nrow(sub)))  # identification of each observation
  n_clusters   <- nlevels(sub$study_id)       # number of studies

  m <- NULL; tipo <- NA
  if (n_clusters >= 2 && nrow(sub) > n_clusters) {
    m <- tryCatch(rma.mv(yi = lnRR, V = v_lnRR, random = ~ 1 | study_id/obs_id,
                         data = sub, method = "REML"),
                  error = function(e) NULL, warning = function(w) NULL)
    if (!is.null(m)) tipo <- "rma.mv"
  }
  if (is.null(m)) {
    m <- tryCatch(rma(yi = lnRR, vi = v_lnRR, data = sub, method = "REML"),
                  error = function(e) NULL)
    if (!is.null(m)) tipo <- "rma"
  }
  if (is.null(m)) return(NULL)

  # Robust variance estimation (RVE): includes dependent effect sizes even when
  # the nature of the dependence is unknown, giving more honest results.
  rob <- NULL
  if (n_clusters >= 3) {
    rob <- tryCatch(robust(m, cluster = sub$study_id), error = function(e) NULL)
  }
  list(model = m, robust = rob, tipo = tipo, n_clusters = n_clusters, data = sub)
}

pull_estimate <- function(ft) {
  src <- if (!is.null(ft$robust)) ft$robust else ft$model
  data.frame(
    Estimate = as.numeric(src$b)[1], SE = as.numeric(src$se)[1],
    Z_val = as.numeric(src$zval)[1], P_val = as.numeric(src$pval)[1],
    CI_lower = as.numeric(src$ci.lb)[1], CI_upper = as.numeric(src$ci.ub)[1],
    robust = !is.null(ft$robust)
  )
}

# ---- omnibus test of a categorical moderator (Q_M) ----
test_cat_moderator <- function(sub, modvar, trait, categoria) {
  vacio <- list(test = NULL, levels = NULL)
  if (!(modvar %in% names(sub))) return(vacio)
  s <- sub[!is.na(sub[[modvar]]), ]
  if (nrow(s) == 0) return(vacio)
  s[[modvar]] <- factor(s[[modvar]])

  est_por_nivel <- tapply(s$study_number, s[[modvar]],
                          function(x) length(unique(x)))
  niveles_ok <- names(est_por_nivel)[!is.na(est_por_nivel) &
                                       est_por_nivel >= MIN_STUDIES_LEVEL]
  if (length(niveles_ok) < MIN_LEVELS) return(vacio)

  s <- s[s[[modvar]] %in% niveles_ok, ]
  s[[modvar]] <- droplevels(factor(s[[modvar]]))
  s$study_id  <- factor(s$study_number)
  s$obs_id    <- factor(seq_len(nrow(s)))

  ajusta <- function(fml) {
    m <- tryCatch(rma.mv(yi = lnRR, V = v_lnRR, mods = fml,
                         random = ~ 1 | study_id/obs_id, data = s, method = "REML"),
                  error = function(e) NULL, warning = function(w) NULL)
    if (is.null(m)) {
      m <- tryCatch(rma(yi = lnRR, vi = v_lnRR, mods = fml, data = s, method = "REML"),
                    error = function(e) NULL)
    }
    m
  }
  m_con <- ajusta(as.formula(paste("~", modvar)))        # WITH intercept -> Q_M test
  m_sin <- ajusta(as.formula(paste("~", modvar, "- 1")))  # NO intercept  -> each level effect
  if (is.null(m_con) || is.null(m_sin)) return(vacio)

  test <- data.frame(
    category = categoria, trait_label = trait, moderator = modvar,
    n_levels = nlevels(s[[modvar]]),
    k_studies = length(unique(s$study_number)), n_obs = m_con$k,
    QM = as.numeric(m_con$QM), QM_df = m_con$m, QM_p = as.numeric(m_con$QMp),
    Signif = signif_stars(as.numeric(m_con$QMp)), stringsAsFactors = FALSE
  )

  nl <- rownames(m_sin$b)
  lev <- sub(paste0("^", modvar), "", nl)
  niveles <- data.frame(
    category = categoria, trait_label = trait, moderator = modvar,
    level = lev, n_obs = as.numeric(table(s[[modvar]])[lev]),
    Estimate = as.numeric(m_sin$b), SE = as.numeric(m_sin$se),
    P_val = as.numeric(m_sin$pval),
    CI_lower = as.numeric(m_sin$ci.lb), CI_upper = as.numeric(m_sin$ci.ub),
    stringsAsFactors = FALSE
  )
  niveles$Signif <- signif_stars(niveles$P_val)
  list(test = test, levels = niveles)
}

# predictor rank: lower = more likely to be the predictor (x) in a pair
pred_rank <- function(t) {
  h <- which(vapply(ASSOC_PREDICTOR_ORDER,
                    function(p) grepl(p, t, ignore.case = TRUE), logical(1)))
  if (length(h)) min(h) else length(ASSOC_PREDICTOR_ORDER) + 1
}
# returns the direction(s) to test for a pair as c(response, predictor)
pair_dirs <- function(a, b) {
  if (isTRUE(ASSOC_BOTH_DIRECTIONS)) return(list(c(a, b), c(b, a)))
  ra <- pred_rank(a); rb <- pred_rank(b)
  if (rb < ra) list(c(a, b)) else list(c(b, a))
}

trait_associations <- function(data_in) {
  d <- data_in[!is.na(data_in$lnRR) & !is.na(data_in$trait_label) &
                 !is.na(data_in$v_lnRR) & data_in$v_lnRR > 0, ]
  if (nrow(d) == 0) return(NULL)
  d$key <- make_unit_key(d)

  traits <- sort(unique(d$trait_label))
  if (length(traits) < 2) return(NULL)

  agg <- aggregate(lnRR ~ key + trait_label, data = d, FUN = mean)
  key_study <- d[!duplicated(d$key), c("key", "study_number")]

  keys <- unique(agg$key)
  M <- matrix(NA_real_, length(keys), length(traits),
              dimnames = list(keys, traits))
  M[cbind(match(agg$key, keys), match(agg$trait_label, traits))] <- agg$lnRR

  cormat <- suppressWarnings(cor(M, use = "pairwise.complete.obs"))

  cor_list <- list()
  for (i in 1:(length(traits) - 1)) {
    for (j in (i + 1):length(traits)) {
      x <- M[, i]; y <- M[, j]
      ok <- is.finite(x) & is.finite(y)
      n_u <- sum(ok)
      if (n_u < 3) next
      n_st <- length(unique(key_study$study_number[key_study$key %in% keys[ok]]))
      if (n_st < ASSOC_MIN_STUDIES) next
      ct <- tryCatch(cor.test(x[ok], y[ok]), error = function(e) NULL)
      if (is.null(ct)) next
      cor_list[[paste(i, j)]] <- data.frame(
        trait_x = traits[i], trait_y = traits[j],
        r = unname(ct$estimate), p = ct$p.value,
        n_units = n_u, n_studies = n_st, stringsAsFactors = FALSE)
    }
  }
  cor_pairs <- bind(cor_list)
  if (nrow(cor_pairs) > 0) {
    cor_pairs$p_adj      <- p.adjust(cor_pairs$p, method = ASSOC_P_ADJUST)
    cor_pairs$Signif     <- signif_stars(cor_pairs$p)
    cor_pairs$Signif_adj <- signif_stars(cor_pairs$p_adj)
    cor_pairs <- cor_pairs[order(-abs(cor_pairs$r)), ]
  }

  agg_split <- split(agg[, c("key", "lnRR")], agg$trait_label)
  mr_list <- list(); pts_list <- list(); fit_list <- list()

  run_pair <- function(resp, pred) {
    ka <- agg_split[[pred]]
    if (is.null(ka)) return(NULL)
    names(ka)[2] <- "pred_x"
    s <- merge(d[d$trait_label == resp, ], ka, by = "key")
    n_st <- length(unique(s$study_number))
    if (nrow(s) < ASSOC_MIN_UNITS || n_st < ASSOC_MIN_STUDIES) return(NULL)
    if (length(unique(s$pred_x)) < 3) return(NULL)
    s$study_id <- factor(s$study_number); s$obs_id <- factor(seq_len(nrow(s)))
    m <- tryCatch(rma.mv(yi = lnRR, V = v_lnRR, mods = ~ pred_x,
                         random = ~ 1 | study_id/obs_id, data = s, method = "REML"),
                  error = function(e) NULL, warning = function(w) NULL)
    tipo <- "rma.mv"
    if (is.null(m)) {
      m <- tryCatch(rma(yi = lnRR, vi = v_lnRR, mods = ~ pred_x,
                        data = s, method = "REML"), error = function(e) NULL)
      tipo <- "rma"
    }
    if (is.null(m) || !("pred_x" %in% rownames(m$b))) return(NULL)

    rob <- if (n_st >= 3) tryCatch(robust(m, cluster = s$study_id),
                                   error = function(e) NULL) else NULL
    src <- if (!is.null(rob)) rob else m

    i <- which(rownames(src$b) == "pred_x")
    g <- seq(min(s$pred_x), max(s$pred_x), length.out = N_GRID)
    pr <- tryCatch(predict(src, newmods = g),
                   error = function(e) predict(m, newmods = g))
    list(
      row = data.frame(trait_y = resp, trait_x = pred, model = tipo,
                       k_studies = n_st, n_obs = m$k,
                       slope = as.numeric(src$b)[i], SE = as.numeric(src$se)[i],
                       P_val = as.numeric(src$pval)[i],
                       CI_lower = src$ci.lb[i], CI_upper = src$ci.ub[i],
                       QM_p = as.numeric(if (!is.null(src$QMp)) src$QMp else m$QMp),
                       R2 = if (tipo == "rma") m$R2 else NA,
                       stringsAsFactors = FALSE),
      pts = data.frame(trait_y = resp, trait_x = pred,
                       x = s$pred_x, y = s$lnRR, weight = 1 / s$v_lnRR,
                       source_file = s$source_file, stringsAsFactors = FALSE),
      fit = data.frame(trait_y = resp, trait_x = pred,
                       x = g, pred = pr$pred, ci_lb = pr$ci.lb, ci_ub = pr$ci.ub,
                       stringsAsFactors = FALSE))
  }

  for (i in 1:(length(traits) - 1)) {
    for (j in (i + 1):length(traits)) {
      for (dir_pair in pair_dirs(traits[i], traits[j])) {
        res <- tryCatch(run_pair(dir_pair[1], dir_pair[2]), error = function(e) NULL)
        if (is.null(res)) next
        tag <- paste(dir_pair[1], "~", dir_pair[2])
        mr_list[[tag]]  <- res$row
        pts_list[[tag]] <- res$pts
        fit_list[[tag]] <- res$fit
      }
    }
  }
  metareg <- bind(mr_list)
  if (nrow(metareg) > 0) {
    metareg$p_adj      <- p.adjust(metareg$P_val, method = ASSOC_P_ADJUST)
    metareg$Signif     <- signif_stars(metareg$P_val)
    metareg$Signif_adj <- signif_stars(metareg$p_adj)
    metareg <- metareg[order(metareg$P_val), ]
  }

  cormat_df <- data.frame(trait = rownames(cormat), cormat,
                          check.names = FALSE, stringsAsFactors = FALSE)

  list(cor_matrix = cormat, cor_matrix_df = cormat_df, cor_pairs = cor_pairs,
       metareg = metareg, points = bind(pts_list), fits = bind(fit_list))
}

ionome_by_genotype <- function(data_in, categoria = IONOME_CATEGORY) {
  if (!("category" %in% names(data_in))) return(NULL)
  cat_ok <- !is.na(data_in$category) &
    tolower(trimws(data_in$category)) == tolower(trimws(categoria))
  d <- data_in[cat_ok & !is.na(data_in[[SPECIE_COL]]) & !is.na(data_in$crop_type) &
                 !is.na(data_in$trait_label) &
                 !is.na(data_in$lnRR) & !is.na(data_in$v_lnRR) & data_in$v_lnRR > 0, ]
  if (nrow(d) == 0) {
    return(NULL)
  }
  if (length(unique(d[[SPECIE_COL]])) < 2) {
    return(NULL)
  }

  minerales <- sort(unique(d$trait_label))
  rows <- list(); glob_rows <- list(); test_list <- list()

  for (tr in minerales) {
    dtr <- d[d$trait_label == tr, ]

    ftg <- tryCatch(fit_trait(dtr), error = function(e) NULL)
    if (!is.null(ftg)) {
      eg <- pull_estimate(ftg)
      glob_rows[[tr]] <- data.frame(
        mineral = tr, crop_type = "(all)", nivel = "GLOBAL (all)",
        k_studies = ftg$n_clusters, n_obs = ftg$model$k,
        lnRR = eg$Estimate, CI_lower = eg$CI_lower, CI_upper = eg$CI_upper,
        pct_change = pct(eg$Estimate), pct_lower = pct(eg$CI_lower),
        pct_upper = pct(eg$CI_upper), P_val = eg$P_val,
        Signif = signif_stars(eg$P_val), stringsAsFactors = FALSE)
    }

    for (cr in sort(unique(dtr$crop_type))) {
      dcr <- dtr[dtr$crop_type == cr, ]
      ref_est <- NA
      ftr <- tryCatch(fit_trait(dcr), error = function(e) NULL)
      if (!is.null(ftr)) ref_est <- pull_estimate(ftr)$Estimate

      for (sp in sort(unique(dcr[[SPECIE_COL]]))) {
        s <- dcr[dcr[[SPECIE_COL]] == sp, ]
        n_st <- length(unique(s$study_number))
        if (nrow(s) < IONOME_MIN_OBS || n_st < IONOME_MIN_STUDIES) next
        ft <- tryCatch(fit_trait(s), error = function(e) NULL)
        if (is.null(ft)) next
        e <- pull_estimate(ft)
        diff0       <- !(e$CI_lower <= 0 & e$CI_upper >= 0)
        mismo_signo <- is.na(ref_est) || sign(e$Estimate) == sign(ref_est)
        rows[[paste(tr, cr, sp)]] <- data.frame(
          mineral = tr, crop_type = cr, nivel = sp,
          k_studies = n_st, n_obs = ft$model$k,
          lnRR = e$Estimate, CI_lower = e$CI_lower, CI_upper = e$CI_upper,
          pct_change = pct(e$Estimate), pct_lower = pct(e$CI_lower),
          pct_upper = pct(e$CI_upper), ref_pct_crop = pct(ref_est),
          P_val = e$P_val, Signif = signif_stars(e$P_val),
          differs_from_zero = diff0,
          follows_trend = (diff0 & mismo_signo) | !diff0,
          stringsAsFactors = FALSE)
      }

      r <- tryCatch(test_cat_moderator(dcr, SPECIE_COL, tr, categoria),
                    error = function(e) NULL)
      if (!is.null(r) && !is.null(r$test)) {
        r$test$crop_type <- cr; test_list[[paste(tr, cr)]] <- r$test
      }
    }
  }

  resumen <- bind(rows)
  if (nrow(resumen) == 0) {
    return(NULL)
  }
  resumen <- resumen[order(resumen$mineral, resumen$crop_type, resumen$pct_change), ]

  resumen$fila <- paste(resumen$crop_type, resumen$nivel, sep = " - ")
  filas <- unique(resumen$fila)
  Mp <- matrix(NA_real_, length(filas), length(minerales),
               dimnames = list(filas, minerales))
  Mp[cbind(match(resumen$fila, filas), match(resumen$mineral, minerales))] <- resumen$pct_change
  matriz <- data.frame(fila = filas, Mp, check.names = FALSE, stringsAsFactors = FALSE)

  combos <- unique(resumen[, c("crop_type", "nivel")])
  gsumm <- list()
  for (i in seq_len(nrow(combos))) {
    cr <- combos$crop_type[i]; sp <- combos$nivel[i]
    x  <- resumen[resumen$crop_type == cr & resumen$nivel == sp, ]
    contra  <- x$mineral[!x$follows_trend]
    dil_sig <- sum(x$differs_from_zero & x$pct_change < 0)
    sub_sig <- sum(x$differs_from_zero & x$pct_change > 0)
    ncontra <- length(contra)
    gsumm[[paste(cr, sp)]] <- data.frame(
      crop_type = cr, genotype = sp,
      n_minerals = nrow(x),
      ionome_mean_pct = round(mean(x$pct_change), 2),
      min_pct = round(min(x$pct_change), 2),
      max_pct = round(max(x$pct_change), 2),
      n_diluted_sig = dil_sig, n_increased_sig = sub_sig,
      n_no_effect = sum(!x$differs_from_zero),
      n_against_trend = ncontra,
      minerals_evaluated   = paste(sort(unique(x$mineral)), collapse = ", "),
      minerals_resisting = paste(sort(contra), collapse = ", "), stringsAsFactors = FALSE)
  }
  resumen_geno <- bind(gsumm)
  if (nrow(resumen_geno) > 0)
    resumen_geno <- resumen_geno[order(-resumen_geno$n_against_trend,
                                       resumen_geno$crop_type,
                                       resumen_geno$ionome_mean_pct), ]

  list(resumen = resumen, resumen_genotipo = resumen_geno, global = bind(glob_rows),
       matriz = matriz, test_por_cultivo = bind(test_list))
}

moderator_summary <- function(data_in) {
  data_in <- data_in[!is.na(data_in$trait_label), ]
  if (nrow(data_in) == 0) return(NULL)
  spt <- tapply(data_in$study_number, data_in$trait_label,
                function(x) length(unique(x)))
  traits <- names(spt)[spt > 1]
  if (length(traits) == 0) return(NULL)

  qm_list <- list(); lv_list <- list()

  one <- function(sub, modvar, trait, categ, crop_lbl) {
    r <- tryCatch(test_cat_moderator(sub, modvar, trait, categ),
                  error = function(e) NULL)
    if (is.null(r) || is.null(r$test)) return(NULL)
    tt <- r$test; tt$crop_type <- crop_lbl
    lv <- NULL
    if (!is.null(r$levels)) {
      lv <- r$levels; lv$crop_type <- crop_lbl
      lv$pct_change <- pct(lv$Estimate)
      lv$pct_lower  <- pct(lv$CI_lower)
      lv$pct_upper  <- pct(lv$CI_upper)
    }
    list(test = tt, levels = lv)
  }

  for (t in traits) {
    sub  <- data_in[data_in$trait_label == t, ]
    categ <- sub$category[1]
    tasks <- list(list(sub, "experiment_type", "(all)"),
                  list(sub, "crop_type",    "(all)"))
    if (SPECIE_COL %in% names(sub))
      for (cr in unique(na.omit(sub$crop_type)))
        tasks[[length(tasks) + 1]] <- list(
          sub[!is.na(sub$crop_type) & sub$crop_type == cr, ], SPECIE_COL, cr)
    for (tk in tasks) {
      res <- one(tk[[1]], tk[[2]], t, categ, tk[[3]])
      if (is.null(res)) next
      qm_list[[length(qm_list) + 1]] <- res$test
      if (!is.null(res$levels)) lv_list[[length(lv_list) + 1]] <- res$levels
    }
  }

  qm <- bind(qm_list); lv <- bind(lv_list)
  if (nrow(qm) == 0) return(NULL)
  qm <- qm[, c("moderator", "category", "trait_label", "crop_type",
               "n_levels", "k_studies", "n_obs", "QM", "QM_df", "QM_p", "Signif")]
  qm <- qm[order(qm$moderator, qm$category, qm$QM_p), ]
  if (nrow(lv) > 0) {
    lv <- lv[, c("moderator", "category", "trait_label", "crop_type", "level",
                 "n_obs", "Estimate", "pct_change", "pct_lower", "pct_upper",
                 "P_val", "Signif")]
    names(lv)[names(lv) == "Estimate"] <- "lnRR"
    lv <- lv[order(lv$moderator, lv$category, lv$trait_label,
                   lv$crop_type, lv$level), ]
  }
  list(qm = qm, levels_paper = lv)
}

repeated_genotype_moderator <- function(data_in) {
  if (!(SPECIE_COL %in% names(data_in))) return(NULL)
  d <- data_in[!is.na(data_in$trait_label) & !is.na(data_in[[SPECIE_COL]]) &
                 !is.na(data_in$crop_type) &
                 !is.na(data_in$lnRR) & !is.na(data_in$v_lnRR) & data_in$v_lnRR > 0, ]
  if (nrow(d) == 0) return(NULL)

  d$geno_id <- paste(d$crop_type, d[[SPECIE_COL]], sep = "|")
  est_por_geno <- tapply(d$study_number, d$geno_id, function(x) length(unique(x)))
  geno_rep <- names(est_por_geno)[est_por_geno >= GENO_MIN_STUDIES]
  if (length(geno_rep) < 2) {
    return(NULL)
  }
  d <- d[d$geno_id %in% geno_rep, ]

  qm_list <- list(); lv_list <- list()

  for (t in unique(d$trait_label)) {
    sub_t <- d[d$trait_label == t, ]; cat_t <- sub_t$category[1]
    for (cr in unique(sub_t$crop_type)) {
      s <- sub_t[sub_t$crop_type == cr, ]
      if (length(unique(s[[SPECIE_COL]])) < 2) next
      r <- tryCatch(test_cat_moderator(s, SPECIE_COL, t, cat_t),
                    error = function(e) NULL)
      if (is.null(r) || is.null(r$test)) next
      r$test$crop_type <- cr
      qm_list[[paste(t, cr)]] <- r$test
      if (!is.null(r$levels)) {
        lv <- r$levels; lv$crop_type <- cr
        lv$pct_change <- pct(lv$Estimate)
        lv$pct_lower  <- pct(lv$CI_lower)
        lv$pct_upper  <- pct(lv$CI_upper)
        lv_list[[paste(t, cr)]] <- lv
      }
    }
  }

  qm <- bind(qm_list); lv <- bind(lv_list)
  if (nrow(qm) == 0) return(NULL)
  qm <- qm[order(qm$category, qm$QM_p), ]
  if (nrow(lv) > 0) {
    lv <- lv[, c("category", "trait_label", "crop_type", "level", "n_obs",
                 "Estimate", "pct_change", "pct_lower", "pct_upper",
                 "P_val", "Signif")]
    names(lv)[names(lv) == "Estimate"] <- "lnRR"
    lv <- lv[order(lv$category, lv$trait_label, lv$crop_type, lv$level), ]
  }
  list(qm = qm, levels = lv, n_genotipos = length(geno_rep),
       genotipos = sort(geno_rep))
}

rawdata <- read_excel(source_file, col_types = "text")
rawdata <- rawdata[, !grepl("^\\.\\.\\.[0-9]+$", names(rawdata)), drop = FALSE]
rawdata <- rawdata[rowSums(!is.na(rawdata) & rawdata != "") > 0, , drop = FALSE]

req_cols <- c("source_file", "trait_label", "category",
              "mean_elevated", "sd_elevated", "nsample",
              "mean_ambient",  "sd_ambient",  "delta_co2_ppm",
              "experiment_type", "crop_type", SPECIE_COL)
faltan <- setdiff(req_cols, names(rawdata))
if (length(faltan) > 0) stop("Missing columns: ", paste(faltan, collapse = ", "))

YEAR_COL <- YEAR_CANDIDATES[YEAR_CANDIDATES %in% names(rawdata)][1]
num_cols <- c("mean_elevated", "sd_elevated", "nsample",
              "mean_ambient", "sd_ambient", "delta_co2_ppm")
if (!is.na(YEAR_COL)) num_cols <- c(num_cols, YEAR_COL)

for (col in num_cols) {
  antes <- sum(is.na(rawdata[[col]]))
  rawdata[[col]] <- to_numeric(rawdata[[col]])
  nuevos <- sum(is.na(rawdata[[col]])) - antes
}

for (col in c("source_file", "trait_label", "category",
              "experiment_type", "crop_type", SPECIE_COL)) {
  rawdata[[col]] <- trimws(rawdata[[col]])
  rawdata[[col]][rawdata[[col]] == ""] <- NA
}

print(sapply(rawdata[num_cols], class))


rawdata$study_number <- as.numeric(as.factor(rawdata$source_file))

# calculate lnRR and v_lnRR
# escalc takes the required data (means, SDs, n sample) and calculates the desired effect size
dat <- escalc(measure = "ROM",
              m1i = mean_elevated, sd1i = sd_elevated, n1i = nsample,
              m2i = mean_ambient,  sd2i = sd_ambient,  n2i = nsample,
              data = rawdata, var.names = c("lnRR", "v_lnRR"))

n_ini <- nrow(dat)
# remove NA effects and negative or zero variance
dat <- dat[!is.na(dat$lnRR) & !is.na(dat$v_lnRR) & dat$v_lnRR > 0, ]
if (nrow(dat) == 0) stop("No valid observations remain.")

run_full_analysis <- function(data_in, out_dir) {

  data_in <- data_in[!is.na(data_in$trait_label), ]
  if (nrow(data_in) == 0) { return(invisible(NULL)) }

  spt <- tapply(data_in$study_number, data_in$trait_label,
                function(x) length(unique(x)))
  data_in <- data_in[data_in$trait_label %in% names(spt)[spt > 1], ]
  if (nrow(data_in) == 0) { return(invisible(NULL)) }

  mods_cat_use <- MODS_CAT
  if (SPECIE_COL %in% mods_cat_use && SPECIE_PARENT %in% names(data_in) &&
      length(unique(na.omit(data_in[[SPECIE_PARENT]]))) > 1) {
    mods_cat_use <- setdiff(mods_cat_use, SPECIE_COL)
  }

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  results_list <- list(); hetero_list <- list(); metareg_list <- list()
  modtest_list <- list(); modlevel_list <- list()

  for (t in unique(data_in$trait_label)) {
    sub   <- subset(data_in, trait_label == t)
    cat_t <- sub$category[1]

    ft <- fit_trait(sub)
    if (is.null(ft)) { next }
    m <- ft$model; est <- pull_estimate(ft)

    results_list[[t]] <- data.frame(
      category = cat_t, trait_label = t, n_obs = m$k, k_studies = ft$n_clusters,
      model = ft$tipo,
      Estimate = est$Estimate, SE = est$SE, Z_val = est$Z_val, P_val = est$P_val,
      CI_lower = est$CI_lower, CI_upper = est$CI_upper,
      pct_change = (exp(est$Estimate) - 1) * 100, stringsAsFactors = FALSE
    )

    if (ft$tipo == "rma.mv") {
      i2 <- i2_multilevel(m, ft$data$v_lnRR)
      hetero_list[[t]] <- data.frame(
        category = cat_t, trait_label = t, model = ft$tipo,
        k_studies = ft$n_clusters, n_obs = m$k,
        sigma2_between = m$sigma2[1],
        sigma2_within = if (length(m$sigma2) > 1) m$sigma2[2] else NA,
        I2_total = i2$I2_total, I2_between = i2$I2_between,
        I2_within = i2$I2_within, Q = m$QE, Q_pval = m$QEp,
        stringsAsFactors = FALSE)
      I2_val <- i2$I2_total
    } else {
      hetero_list[[t]] <- data.frame(
        category = cat_t, trait_label = t, model = ft$tipo,
        k_studies = ft$n_clusters, n_obs = m$k,
        sigma2_between = m$tau2, sigma2_within = NA,
        I2_total = m$I2, I2_between = NA, I2_within = NA,
        Q = m$QE, Q_pval = m$QEp, stringsAsFactors = FALSE)
      I2_val <- m$I2
    }

    for (mv in mods_cat_use) {
      r <- test_cat_moderator(sub, mv, t, cat_t)
      if (!is.null(r$test))   modtest_list[[paste(t, mv)]]  <- r$test
      if (!is.null(r$levels)) modlevel_list[[paste(t, mv)]] <- r$levels
    }

    sub_ok <- sub[!is.na(sub$delta_co2_ppm), ]
    vacia <- data.frame(category = cat_t, trait_label = t,
                        k_studies = ft$n_clusters, n_obs = m$k, I2 = I2_val,
                        adjusted = FALSE, Estimate = NA, SE = NA, Z_val = NA,
                        P_val = NA, CI_lower = NA, CI_upper = NA,
                        QM = NA, QM_p = NA, R2 = NA, stringsAsFactors = FALSE)
    if (nrow(sub_ok) < 4 || length(unique(sub_ok$delta_co2_ppm)) < 3) {
      metareg_list[[t]] <- vacia; next
    }

    sub_ok$study_id <- factor(sub_ok$study_number)
    sub_ok$obs_id   <- factor(seq_len(nrow(sub_ok)))
    # meta-regression: does the CO2 effect on a trait depend on the DOSE (delta CO2)?
    mr <- tryCatch(rma.mv(yi = lnRR, V = v_lnRR, mods = ~ delta_co2_ppm,   # continuous predictor (the slope)
                          random = ~ 1 | study_id/obs_id, data = sub_ok, method = "REML"),
                   error = function(e) NULL, warning = function(w) NULL)
    tipo_mr <- "rma.mv"
    if (is.null(mr)) {
      mr <- tryCatch(rma(yi = lnRR, vi = v_lnRR, mods = ~ delta_co2_ppm,
                         data = sub_ok, method = "REML"), error = function(e) NULL)
      tipo_mr <- "rma"
    }
    if (is.null(mr) || !("delta_co2_ppm" %in% rownames(mr$b))) {
      metareg_list[[t]] <- vacia; next
    }

    idx <- which(rownames(mr$b) == "delta_co2_ppm")
    R2_val <- if (tipo_mr == "rma") mr$R2 else NA
    if (tipo_mr == "rma.mv") {
      # pseudo-R2: NULL model (no predictor) vs model with predictor
      m0 <- tryCatch(rma.mv(yi = lnRR, V = v_lnRR, random = ~ 1 | study_id/obs_id,
                            data = sub_ok, method = "REML"),
                     error = function(e) NULL, warning = function(w) NULL)
      if (!is.null(m0)) {
        s0 <- sum(m0$sigma2); s1 <- sum(mr$sigma2)
        if (is.finite(s0) && s0 > 0) R2_val <- max(0, 100 * (s0 - s1) / s0)
      }
    }
    metareg_list[[t]] <- data.frame(
      category = cat_t, trait_label = t,
      k_studies = length(unique(sub_ok$study_number)), n_obs = mr$k,
      I2 = I2_val, adjusted = TRUE,
      Estimate = as.numeric(mr$b)[idx], SE = as.numeric(mr$se)[idx],
      Z_val = as.numeric(mr$zval)[idx], P_val = as.numeric(mr$pval)[idx],
      CI_lower = mr$ci.lb[idx], CI_upper = mr$ci.ub[idx],
      QM = as.numeric(mr$QM), QM_p = as.numeric(mr$QMp),
      R2 = R2_val, stringsAsFactors = FALSE)

  }

  if (length(results_list) == 0) { return(invisible(NULL)) }

  results_model <- bind(results_list)
  results_model <- results_model[order(results_model$category, results_model$trait_label), ]
  results_model$Signif <- signif_stars(results_model$P_val)
  heterogenity <- bind(hetero_list); metareg_model <- bind(metareg_list)
  mod_tests <- bind(modtest_list);   mod_levels <- bind(modlevel_list)

  summary_by_category <- results_model[, c(
    "category", "trait_label", "k_studies", "n_obs", "model",
    "Estimate", "SE", "CI_lower", "CI_upper", "pct_change",
    "Z_val", "P_val", "Signif")]

  assoc <- if (MAKE_ASSOC) tryCatch(trait_associations(data_in),
                                    error = function(e) NULL) else NULL
  iono  <- if (MAKE_IONOME_BY_GENO)
    tryCatch(ionome_by_genotype(data_in),
             error = function(e) NULL)
  else NULL
  modsum <- tryCatch(moderator_summary(data_in),
                     error = function(e) NULL)
  genorep <- if (MAKE_REPEATED_GENO)
    tryCatch(repeated_genotype_moderator(data_in),
             error = function(e) NULL)
  else NULL

  hojas <- list("Results" = results_model, "Summary_by_category" = summary_by_category,
                "Heterogenity" = heterogenity, "MetaRegression" = metareg_model,
                "Moderator_tests" = mod_tests, "Moderator_levels" = mod_levels)
  if (!is.null(assoc)) {
    if (nrow(assoc$cor_pairs) > 0) hojas[["Trait_correlations"]]    <- assoc$cor_pairs
    if (nrow(assoc$metareg)   > 0) hojas[["Trait_metaregressions"]] <- assoc$metareg
  }
  if (!is.null(iono)) {
    if (!is.null(iono$resumen_genotipo) && nrow(iono$resumen_genotipo) > 0)
      hojas[["Ionome_genotype_summary"]] <- iono$resumen_genotipo
    hojas[["Ionome_by_genotype"]] <- iono$resumen
    if (!is.null(iono$test_por_cultivo) && nrow(iono$test_por_cultivo) > 0)
      hojas[["Ionome_geno_test"]] <- iono$test_por_cultivo
  }
  if (!is.null(modsum)) {
    hojas[["Moderator_summary"]] <- modsum$qm
    if (!is.null(modsum$levels_paper) && nrow(modsum$levels_paper) > 0)
      hojas[["Moderator_levels_paper"]] <- modsum$levels_paper
  }
  if (!is.null(genorep)) {
    hojas[["RepGeno_moderator"]] <- genorep$qm
    if (!is.null(genorep$levels) && nrow(genorep$levels) > 0)
      hojas[["RepGeno_levels"]] <- genorep$levels
  }
  write_xlsx(lapply(hojas, function(x) if (is.null(x)) data.frame() else x),
             path = file.path(out_dir, "results.xlsx"))

  invisible(NULL)
}

run_full_analysis(dat, file.path(base_dir, "overall"))

for (et in unique(na.omit(dat$experiment_type))) {
  run_full_analysis(dat[!is.na(dat$experiment_type) & dat$experiment_type == et, ],
                    file.path(base_dir, "by_experiment_type", clean_name(et)))
}

for (cs in unique(na.omit(dat$crop_type))) {
  run_full_analysis(dat[!is.na(dat$crop_type) & dat$crop_type == cs, ],
                    file.path(base_dir, "by_crop_type", clean_name(cs)))
}

for (sp in unique(na.omit(dat[[SPECIE_COL]]))) {
  run_full_analysis(dat[!is.na(dat[[SPECIE_COL]]) & dat[[SPECIE_COL]] == sp, ],
                    file.path(base_dir, "by_genotype", clean_name(sp)))
}

combos <- unique(dat[, c("experiment_type", "crop_type")])
combos <- combos[!is.na(combos$experiment_type) & !is.na(combos$crop_type), ]
for (i in seq_len(nrow(combos))) {
  et <- combos$experiment_type[i]; cs <- combos$crop_type[i]
  run_full_analysis(
    dat[!is.na(dat$experiment_type) & !is.na(dat$crop_type) &
          dat$experiment_type == et & dat$crop_type == cs, ],
    file.path(base_dir, "by_experiment_type_crop_species",
              paste0(clean_name(et), "__", clean_name(cs))))
}

combos_cs <- unique(dat[, c("crop_type", SPECIE_COL)])
combos_cs <- combos_cs[!is.na(combos_cs$crop_type) & !is.na(combos_cs[[SPECIE_COL]]), ]
for (i in seq_len(nrow(combos_cs))) {
  cs <- combos_cs$crop_type[i]; sp <- combos_cs[[SPECIE_COL]][i]
  run_full_analysis(
    dat[!is.na(dat$crop_type) & !is.na(dat[[SPECIE_COL]]) &
          dat$crop_type == cs & dat[[SPECIE_COL]] == sp, ],
    file.path(base_dir, "by_crop_type_genotype",
              paste0(clean_name(cs), "__", clean_name(sp))))
}
