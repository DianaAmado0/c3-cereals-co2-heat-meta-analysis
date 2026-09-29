library(readxl)
library(metafor)
library(writexl)


if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}


source_file <- "hs_dataset.xlsx"
out_dir    <- "RESULTS/HS"

CO2_COL   <- "ambient_co2_ppm" 
TEMP_COL  <- "t_elevated"       
SPECIE_COL    <- "genotype"
SPECIE_PARENT <- "crop_type" 

MODS_CAT <- c("co2_level", "experiment_type", "crop_type", SPECIE_COL)
MIN_STUDIES_LEVEL <- 2
MIN_LEVELS        <- 2

DIGITS <- 6; N_GRID <- 100

clean_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))

signif_stars <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*",
                              ifelse(p < 0.1, ".", "ns")))))
}

pct <- function(x) (exp(x) - 1) * 100

to_numeric <- function(x) {
  if (is.numeric(x)) return(x)
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "na", "n/a", "n/d", "-", "--", ".")] <- NA
  x <- gsub("\\s", "", x); x <- gsub(",", ".", x, fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

bind <- function(l) if (length(l) > 0) do.call(rbind, l) else data.frame()

# ---- heterogeneity (I2), partitioned into between/within studies ----
i2_multilevel <- function(m, vi) {
  out <- list(I2_total = NA, I2_between = NA, I2_within = NA)
  res <- tryCatch({
    W <- diag(1 / vi); X <- model.matrix(m)
    P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
    s2 <- (m$k - m$p) / sum(diag(P)); tot <- sum(m$sigma2)
    list(I2_total = 100 * tot / (tot + s2),
         I2_between = 100 * m$sigma2[1] / (tot + s2),
         I2_within = if (length(m$sigma2) > 1) 100 * m$sigma2[2] / (tot + s2) else NA)
  }, error = function(e) NULL)
  if (!is.null(res)) out <- res
  out
}

# ---- multilevel model fit (applied per trait) ----
fit_trait <- function(sub) {
  sub$study_id <- factor(sub$study_number); sub$obs_id <- factor(seq_len(nrow(sub)))
  n_clusters <- nlevels(sub$study_id); m <- NULL; tipo <- NA
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
  rob <- if (n_clusters >= 3) tryCatch(robust(m, cluster = sub$study_id),
                                       error = function(e) NULL) else NULL
  list(model = m, robust = rob, tipo = tipo, n_clusters = n_clusters, data = sub)
}

pull_estimate <- function(ft) {
  src <- if (!is.null(ft$robust)) ft$robust else ft$model
  data.frame(Estimate = as.numeric(src$b)[1], SE = as.numeric(src$se)[1],
             P_val = as.numeric(src$pval)[1], CI_lower = as.numeric(src$ci.lb)[1],
             CI_upper = as.numeric(src$ci.ub)[1], robust = !is.null(ft$robust))
}

# ---- omnibus test of a categorical moderator (Q_M) ----
test_cat_moderator <- function(sub, modvar, trait, categoria) {
  vacio <- list(test = NULL, levels = NULL)
  if (!(modvar %in% names(sub))) return(vacio)
  s <- sub[!is.na(sub[[modvar]]), ]; if (nrow(s) == 0) return(vacio)
  s[[modvar]] <- factor(s[[modvar]])
  est_x_niv <- tapply(s$study_number, s[[modvar]], function(x) length(unique(x)))
  ok_niv <- names(est_x_niv)[!is.na(est_x_niv) & est_x_niv >= MIN_STUDIES_LEVEL]
  if (length(ok_niv) < MIN_LEVELS) return(vacio)
  s <- s[s[[modvar]] %in% ok_niv, ]; s[[modvar]] <- droplevels(factor(s[[modvar]]))
  s$study_id <- factor(s$study_number); s$obs_id <- factor(seq_len(nrow(s)))
  ajusta <- function(fml) {
    m <- tryCatch(rma.mv(yi = lnRR, V = v_lnRR, mods = fml,
                         random = ~ 1 | study_id/obs_id, data = s, method = "REML"),
                  error = function(e) NULL, warning = function(w) NULL)
    if (is.null(m)) m <- tryCatch(rma(yi = lnRR, vi = v_lnRR, mods = fml,
                                      data = s, method = "REML"), error = function(e) NULL)
    m
  }
  m_con <- ajusta(as.formula(paste("~", modvar)))
  m_sin <- ajusta(as.formula(paste("~", modvar, "- 1")))
  if (is.null(m_con) || is.null(m_sin)) return(vacio)
  test <- data.frame(category = categoria, trait_label = trait, moderator = modvar,
                     n_levels = nlevels(s[[modvar]]),
                     k_studies = length(unique(s$study_number)), n_obs = m_con$k,
                     QM = as.numeric(m_con$QM), QM_df = m_con$m, QM_p = as.numeric(m_con$QMp),
                     Signif = signif_stars(as.numeric(m_con$QMp)), stringsAsFactors = FALSE)
  lev <- sub(paste0("^", modvar), "", rownames(m_sin$b))
  niveles <- data.frame(category = categoria, trait_label = trait, moderator = modvar,
                        level = lev, n_obs = as.numeric(table(s[[modvar]])[lev]),
                        Estimate = as.numeric(m_sin$b),
                        pct_change = pct(as.numeric(m_sin$b)),
                        pct_lower = pct(as.numeric(m_sin$ci.lb)),
                        pct_upper = pct(as.numeric(m_sin$ci.ub)),
                        P_val = as.numeric(m_sin$pval),
                        Signif = signif_stars(as.numeric(m_sin$pval)),
                        stringsAsFactors = FALSE)
  list(test = test, levels = niveles)
}

raw <- read_excel(source_file, col_types = "text")
raw <- raw[, !grepl("^\\.\\.\\.[0-9]+$", names(raw)), drop = FALSE]
raw <- raw[rowSums(!is.na(raw) & raw != "") > 0, , drop = FALSE]

if ("n_sample" %in% names(raw) && !("nsample" %in% names(raw)))
  names(raw)[names(raw) == "n_sample"] <- "nsample"

req <- c("trait_label", "category", "mean_ambient", "sd_ambient",
         "mean_elevated", "sd_elevated", "nsample",
         "experiment_type", "crop_type", SPECIE_COL, CO2_COL, TEMP_COL)
faltan <- setdiff(req, names(raw))
if (length(faltan) > 0) stop("Missing columns: ", paste(faltan, collapse = ", "))

num_cols <- c("mean_ambient", "sd_ambient", "mean_elevated", "sd_elevated",
              "nsample", CO2_COL, TEMP_COL)
for (col in num_cols) raw[[col]] <- to_numeric(raw[[col]])

for (col in c("trait_label", "category", "experiment_type", "crop_type", SPECIE_COL)) {
  raw[[col]] <- trimws(raw[[col]]); raw[[col]][raw[[col]] == ""] <- NA
}

raw$co2_level <- ifelse(is.na(raw[[CO2_COL]]), NA, paste0(round(raw[[CO2_COL]]), "ppm"))

raw$study_number <- as.numeric(as.factor(raw$source_file))

# calculate lnRR and v_lnRR (heat effect: elevated temp vs ambient temp)
dat <- escalc(measure = "ROM",
              m1i = mean_elevated, sd1i = sd_elevated, n1i = nsample,
              m2i = mean_ambient,  sd2i = sd_ambient,  n2i = nsample,
              data = raw, var.names = c("lnRR", "v_lnRR"))
n_ini <- nrow(dat)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

spt <- tapply(dat$study_number, dat$trait_label, function(x) length(unique(x)))
dat <- dat[dat$trait_label %in% names(spt)[spt > 1], ]

res_list <- list(); het_list <- list(); mr_list <- list()

for (t in unique(dat$trait_label)) {
  sub <- dat[dat$trait_label == t, ]; cat_t <- sub$category[1]
  ft <- fit_trait(sub); if (is.null(ft)) next
  m <- ft$model; e <- pull_estimate(ft)
  
  res_list[[t]] <- data.frame(
    category = cat_t, trait_label = t, k_studies = ft$n_clusters, n_obs = m$k,
    modelo = ft$tipo,
    Estimate = e$Estimate, CI_lower = e$CI_lower, CI_upper = e$CI_upper,
    pct_change = pct(e$Estimate), pct_lower = pct(e$CI_lower), pct_upper = pct(e$CI_upper),
    P_val = e$P_val, Signif = signif_stars(e$P_val), stringsAsFactors = FALSE)
  
  if (ft$tipo == "rma.mv") {
    i2 <- i2_multilevel(m, ft$data$v_lnRR)
    het_list[[t]] <- data.frame(category = cat_t, trait_label = t,
                                k_studies = ft$n_clusters, n_obs = m$k, I2_total = i2$I2_total,
                                I2_between = i2$I2_between, I2_within = i2$I2_within,
                                Q = m$QE, Q_pval = m$QEp, stringsAsFactors = FALSE)
  } else {
    het_list[[t]] <- data.frame(category = cat_t, trait_label = t,
                                k_studies = ft$n_clusters, n_obs = m$k, I2_total = m$I2,
                                I2_between = NA, I2_within = NA, Q = m$QE, Q_pval = m$QEp, stringsAsFactors = FALSE)
  }
  
  # meta-regression vs t_elevated (heat magnitude)
  s2 <- sub[!is.na(sub[[TEMP_COL]]), ]
  if (nrow(s2) >= 4 && length(unique(s2[[TEMP_COL]])) >= 3) {
    s2$tmp <- s2[[TEMP_COL]]; s2$study_id <- factor(s2$study_number)
    s2$obs_id <- factor(seq_len(nrow(s2)))
    mr <- tryCatch(rma.mv(yi = lnRR, V = v_lnRR, mods = ~ tmp,
                          random = ~ 1 | study_id/obs_id, data = s2, method = "REML"),
                   error = function(e) NULL, warning = function(w) NULL)
    tipo_mr <- "rma.mv"
    if (is.null(mr)) { mr <- tryCatch(rma(yi = lnRR, vi = v_lnRR, mods = ~ tmp,
                                          data = s2, method = "REML"), error = function(e) NULL)
    tipo_mr <- "rma" }
    if (!is.null(mr) && "tmp" %in% rownames(mr$b)) {
      rob <- if (length(unique(s2$study_number)) >= 3)
        tryCatch(robust(mr, cluster = s2$study_id), error = function(e) NULL) else NULL
      src <- if (!is.null(rob)) rob else mr
      i <- which(rownames(src$b) == "tmp")
      # pseudo-R2
      R2v <- if (tipo_mr == "rma") mr$R2 else NA
      if (tipo_mr == "rma.mv") {
        m0 <- tryCatch(rma.mv(yi = lnRR, V = v_lnRR, random = ~ 1 | study_id/obs_id,
                              data = s2, method = "REML"), error = function(e) NULL,
                       warning = function(w) NULL)
        if (!is.null(m0)) { s0 <- sum(m0$sigma2); s1 <- sum(mr$sigma2)
        if (is.finite(s0) && s0 > 0) R2v <- max(0, 100 * (s0 - s1) / s0) }
      }
      mr_list[[t]] <- data.frame(category = cat_t, trait_label = t,
                                 k_studies = length(unique(s2$study_number)), n_obs = mr$k,
                                 slope = as.numeric(src$b)[i], SE = as.numeric(src$se)[i],
                                 P_val = as.numeric(src$pval)[i], CI_lower = src$ci.lb[i], CI_upper = src$ci.ub[i],
                                 R2_pct = R2v, Signif = signif_stars(as.numeric(src$pval)[i]), stringsAsFactors = FALSE)
    }
  }
}

results   <- bind(res_list); results <- results[order(results$category, results$trait_label), ]
hetero    <- bind(het_list)
metareg_t <- bind(mr_list)

qm_list <- list(); lv_list <- list()
one <- function(s, modvar, trait, categ, crop_lbl) {
  r <- tryCatch(test_cat_moderator(s, modvar, trait, categ), error = function(e) NULL)
  if (is.null(r) || is.null(r$test)) return(NULL)
  tt <- r$test; tt$crop_type <- crop_lbl
  lv <- r$levels; if (!is.null(lv)) lv$crop_type <- crop_lbl
  list(test = tt, levels = lv)
}
for (t in unique(dat$trait_label)) {
  sub <- dat[dat$trait_label == t, ]; categ <- sub$category[1]
  tasks <- list(list(sub, "co2_level", "(all)"),
                list(sub, "experiment_type", "(all)"),
                list(sub, "crop_type", "(all)"))
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
mod_qm <- bind(qm_list); mod_lv <- bind(lv_list)
if (nrow(mod_qm) > 0) mod_qm <- mod_qm[order(mod_qm$moderator, mod_qm$category, mod_qm$QM_p), ]

mitig_list <- list(); split_list <- list()

co2_num <- suppressWarnings(as.numeric(gsub("ppm", "", dat$co2_level)))
dat$co2num <- co2_num
dat$co2cat <- ave(dat$co2num, dat$study_number, FUN = function(x) {
  if (length(unique(x[!is.na(x)])) < 2) return(rep(NA_character_, length(x)))
  ifelse(x <= min(x, na.rm = TRUE), "aCO2",
         ifelse(x >= max(x, na.rm = TRUE), "eCO2", NA_character_))
})
if (sum(!is.na(dat$co2cat)) > 0) {
  dat$unit_id <- paste(dat$study_number, dat$crop_type, dat[[SPECIE_COL]], sep = "|")
  
  for (t in unique(dat$trait_label)) {
    sub <- dat[dat$trait_label == t, ]; cat_t <- sub$category[1]
    
    for (cc in c("aCO2", "eCO2")) {
      sc <- sub[!is.na(sub$co2cat) & sub$co2cat == cc, ]
      if (nrow(sc) == 0) next
      ftc <- tryCatch(fit_trait(sc), error = function(e) NULL); if (is.null(ftc)) next
      ec <- pull_estimate(ftc)
      split_list[[paste(t, cc)]] <- data.frame(
        category = cat_t, trait_label = t,
        co2_level = ifelse(cc == "aCO2", "ambient CO2", "elevated CO2"),
        k_studies = ftc$n_clusters, n_obs = ftc$model$k,
        pct_change = pct(ec$Estimate), pct_lower = pct(ec$CI_lower),
        pct_upper = pct(ec$CI_upper), P_val = ec$P_val,
        Signif = signif_stars(ec$P_val), stringsAsFactors = FALSE)
    }
    
    agg <- aggregate(lnRR ~ unit_id + co2cat, data = sub, FUN = mean)
    wide <- reshape(agg, idvar = "unit_id", timevar = "co2cat", direction = "wide")
    col_lo <- "lnRR.aCO2"; col_hi <- "lnRR.eCO2"
    if (!all(c(col_lo, col_hi) %in% names(wide))) next
    wide$diff <- wide[[col_hi]] - wide[[col_lo]]
    wide <- wide[is.finite(wide$diff), ]
    n_pares <- nrow(wide)
    n_est <- length(unique(sub$study_number[sub$unit_id %in% wide$unit_id]))
    if (n_pares < 3) next
    
    vsum <- aggregate(v_lnRR ~ unit_id, data = sub, FUN = function(x) sum(x) / length(x))
    wide <- merge(wide, vsum, by = "unit_id")
    wide$study_id <- factor(sub$study_number[match(wide$unit_id, sub$unit_id)])
    wide$obs_id <- factor(seq_len(nrow(wide)))
    
    md <- tryCatch(rma.mv(yi = diff, V = v_lnRR, random = ~ 1 | study_id/obs_id,
                          data = wide, method = "REML"),
                   error = function(e) NULL, warning = function(w) NULL)
    if (is.null(md)) md <- tryCatch(rma(yi = diff, vi = v_lnRR, data = wide, method = "REML"),
                                    error = function(e) NULL)
    if (is.null(md)) next
    rob <- if (n_est >= 3) tryCatch(robust(md, cluster = wide$study_id),
                                    error = function(e) NULL) else NULL
    src <- if (!is.null(rob)) rob else md
    est <- as.numeric(src$b)[1]; p <- as.numeric(src$pval)[1]
    lb <- src$ci.lb[1]; ub <- src$ci.ub[1]
    
    mitig_list[[t]] <- data.frame(
      category = cat_t, trait_label = t, n_pairs = n_pares, k_studies = n_est,
      diff_lnRR = est, diff_pct = pct(est), CI_lower_pct = pct(lb), CI_upper_pct = pct(ub),
      P_val = p, Signif = signif_stars(p), stringsAsFactors = FALSE)
  }
}
mitig <- bind(mitig_list); heat_split <- bind(split_list)
if (nrow(mitig) > 0) mitig <- mitig[order(mitig$P_val), ]

hojas <- list(Results = results, Heterogenity = hetero, MetaRegression_temp = metareg_t,
              Moderator_summary = mod_qm, Moderator_levels_paper = mod_lv,
              CO2_mitigation = mitig, Heat_effect_by_CO2 = heat_split)
write_xlsx(lapply(hojas, function(x) if (is.null(x) || nrow(x) == 0) data.frame() else x),
           path = file.path(out_dir, "results.xlsx"))