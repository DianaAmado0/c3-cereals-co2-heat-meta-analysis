
# SENSITIVITY ANALYSIS - risk of bias
# Full results are compared with results after removing high-risk studies.

library(readxl)
library(writexl)

if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

full_path <- "RESULTS/overall/results.xlsx"
bias_path <- "RESULTS BIAS/overall/results.xlsx"
out_path  <- "RESULTS/sensitivity_bias_comparison.xlsx"

full <- as.data.frame(read_excel(full_path, sheet = "Results"))
bias <- as.data.frame(read_excel(bias_path, sheet = "Results"))

keep <- c("category", "trait_label", "k_studies",
          "Estimate", "CI_lower", "CI_upper", "P_val")

miss_full <- setdiff(keep, names(full))
miss_bias <- setdiff(keep, names(bias))
if (length(miss_full)) stop("Full results missing columns: ", paste(miss_full, collapse = ", "))
if (length(miss_bias)) stop("BIAS results missing columns: ", paste(miss_bias, collapse = ", "))

pct <- function(x) (exp(x) - 1) * 100

comp <- merge(full[keep], bias[keep],
              by = c("category", "trait_label"),
              suffixes = c("_all", "_bias"))

comp$pct_all        <- round(pct(comp$Estimate_all), 2)
comp$pct_bias       <- round(pct(comp$Estimate_bias), 2)
comp$delta_lnRR     <- round(comp$Estimate_bias - comp$Estimate_all, 4)
comp$pct_change_rel <- round(100 * (comp$Estimate_bias - comp$Estimate_all) / abs(comp$Estimate_all), 1)
comp$k_removed      <- comp$k_studies_all - comp$k_studies_bias

comp <- comp[order(-abs(comp$pct_change_rel)), ]

# keep only value columns
comp <- comp[, c("category", "trait_label",
                 "k_studies_all", "k_studies_bias", "k_removed",
                 "Estimate_all", "P_val_all", "pct_all",
                 "Estimate_bias", "P_val_bias", "pct_bias",
                 "delta_lnRR", "pct_change_rel")]

print(comp, row.names = FALSE)

cats <- sort(unique(comp$category))
cat_summary <- do.call(rbind, lapply(cats, function(ct) {
  d <- comp[comp$category == ct, ]
  data.frame(
    category            = ct,
    n_traits            = nrow(d),
    studies_removed_max = max(d$k_removed, na.rm = TRUE),
    median_abs_change   = round(median(abs(d$pct_change_rel), na.rm = TRUE), 1),
    max_abs_change      = round(max(abs(d$pct_change_rel),    na.rm = TRUE), 1),
    stringsAsFactors = FALSE
  )
}))

print(cat_summary, row.names = FALSE)

write_xlsx(
  list(trait_level = comp, category_summary = cat_summary),
  path = out_path
)
cat("\nSaved to:", normalizePath(out_path), "\n")