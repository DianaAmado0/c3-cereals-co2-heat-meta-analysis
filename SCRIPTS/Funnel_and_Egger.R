# FUNNEL PLOTS and EGGER'S TEST (small-study effects) for publication BIAS

library(readxl)
library(metafor)
library(ggplot2)
library(dplyr)
library(writexl)

if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

source_file  <- "dataset.xlsx"
results_xlsx <- "RESULTS/overall/results.xlsx"
output_dir   <- "RESULTS/overall/funnel"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

MIN_STUDIES_EGGER <- 4

to_numeric <- function(x) {
  if (is.numeric(x)) return(x)
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "na", "n/a", "n/d", "-", "--", ".")] <- NA
  x <- gsub("\\s", "", x); x <- gsub(",", ".", x, fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

signif_stars <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*",
                              ifelse(p < 0.1, ".", "ns")))))
}

raw <- read_excel(source_file, col_types = "text")
raw <- raw[, !grepl("^\\.\\.\\.[0-9]+$", names(raw)), drop = FALSE]
raw <- raw[rowSums(!is.na(raw) & raw != "") > 0, , drop = FALSE]
for (col in c("mean_elevated", "sd_elevated", "nsample", "mean_ambient", "sd_ambient"))
  raw[[col]] <- to_numeric(raw[[col]])
raw$trait_label <- trimws(raw$trait_label)
raw$category    <- trimws(raw$category)
raw$study_number <- as.numeric(as.factor(raw$source_file))

dat_es <- escalc(measure = "ROM",
                 m1i = mean_elevated, sd1i = sd_elevated, n1i = nsample,
                 m2i = mean_ambient,  sd2i = sd_ambient,  n2i = nsample,
                 data = raw, var.names = c("lnRR", "v_lnRR"))
dat_es <- dat_es[!is.na(dat_es$lnRR) & !is.na(dat_es$v_lnRR) & dat_es$v_lnRR > 0, ]
dat_es$yi <- dat_es$lnRR
dat_es$SE <- sqrt(dat_es$v_lnRR)


sheets <- excel_sheets(results_xlsx)
res_sheet <- if ("Results" %in% sheets) "Results" else sheets[1]
results <- as.data.frame(read_excel(results_xlsx, sheet = res_sheet))

results <- results[results$trait_label %in% unique(dat_es$trait_label), ]

panel_info <- results %>%
  transmute(category = category,
            trait_label = trait_label,
            k_studies = k_studies,
            pooled_lnRR = Estimate) %>%
  arrange(category, trait_label)

trait_order <- panel_info$trait_label
dat_es$trait_label     <- factor(dat_es$trait_label, levels = trait_order)
panel_info$trait_label <- factor(panel_info$trait_label, levels = trait_order)
dat_es <- dat_es[!is.na(dat_es$trait_label), ]


egger_rows <- list()
for (t in trait_order) {
  s <- dat_es[dat_es$trait_label == t, ]
  n_st <- length(unique(s$study_number))
  if (n_st < MIN_STUDIES_EGGER) next
  s$study_id <- factor(s$study_number); s$obs_id <- factor(seq_len(nrow(s)))
  s$sei <- sqrt(s$v_lnRR)
  eg <- tryCatch(rma.mv(yi = lnRR, V = v_lnRR, mods = ~ sei,
                        random = ~ 1 | study_id/obs_id, data = s, method = "REML"),
                 error = function(e) NULL, warning = function(w) NULL)
  if (is.null(eg))
    eg <- tryCatch(rma(yi = lnRR, vi = v_lnRR, mods = ~ sei, data = s, method = "REML"),
                   error = function(e) NULL)
  if (is.null(eg) || !("sei" %in% rownames(eg$b))) next
  i <- which(rownames(eg$b) == "sei")
  egger_rows[[t]] <- data.frame(
    category    = as.character(panel_info$category[panel_info$trait_label == t][1]),
    trait_label = t,
    k_studies   = n_st,
    n_obs       = eg$k,
    egger_slope = as.numeric(eg$b)[i],
    egger_SE    = as.numeric(eg$se)[i],
    egger_p     = as.numeric(eg$pval)[i],
    Signif      = signif_stars(as.numeric(eg$pval)[i]),
    stringsAsFactors = FALSE)
}
egger_table <- if (length(egger_rows) > 0) do.call(rbind, egger_rows) else data.frame()
if (nrow(egger_table) > 0) egger_table <- egger_table[order(egger_table$egger_p), ]

write_xlsx(list(Egger_test = egger_table),
           path = file.path(output_dir, "egger_test.xlsx"))
write.csv(egger_table, file.path(output_dir, "egger_test.csv"), row.names = FALSE)


panel_info$panel_letter <- paste0(letters[seq_len(nrow(panel_info))], ")")
panel_info$facet_title <- paste0(panel_info$panel_letter, "  ",
                                 as.character(panel_info$trait_label),
                                 "\n", "k = ", panel_info$k_studies)
facet_names <- setNames(panel_info$facet_title, as.character(panel_info$trait_label))

reference_lines <- panel_info %>%
  filter(!is.na(pooled_lnRR)) %>%
  select(trait_label, pooled_lnRR)

p_funnel <- ggplot(dat_es, aes(x = yi, y = SE)) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.45, colour = "black") +
  geom_vline(data = reference_lines, aes(xintercept = pooled_lnRR),
             inherit.aes = FALSE, linetype = "solid", linewidth = 0.65, colour = "black") +
  geom_point(shape = 21, size = 2.3, stroke = 0.60, fill = "white",
             colour = "black", alpha = 0.90) +
  scale_y_reverse(expand = expansion(mult = c(0.08, 0.08))) +
  facet_wrap(~ trait_label, ncol = 4, scales = "free",
             labeller = as_labeller(facet_names)) +
  labs(x = "Log response ratio (lnRR)", y = "Standard error") +
  theme_classic(base_family = "Arial", base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(family = "Arial", face = "bold", size = 10,
                              colour = "black", hjust = 0, lineheight = 1.10,
                              margin = margin(t = 4, b = 8)),
    axis.title = element_text(family = "Arial", face = "bold", size = 12, colour = "black"),
    axis.text  = element_text(family = "Arial", size = 9, colour = "black"),
    axis.line  = element_line(linewidth = 0.55, colour = "black"),
    axis.ticks = element_line(linewidth = 0.45, colour = "black"),
    axis.ticks.length = unit(1.8, "mm"),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.spacing.x = unit(2.0, "lines"),
    panel.spacing.y = unit(2.2, "lines"),
    plot.margin = margin(15, 15, 15, 15)
  )

n_panels  <- nrow(panel_info)
n_columns <- 4
n_rows    <- ceiling(n_panels / n_columns)
figure_width  <- 15
figure_height <- n_rows * 3.6


ggsave(file.path(output_dir, "FigureFunnel.tiff"),
       p_funnel, width = figure_width, height = figure_height,
       units = "in", dpi = 600, compression = "lzw", limitsize = FALSE)

