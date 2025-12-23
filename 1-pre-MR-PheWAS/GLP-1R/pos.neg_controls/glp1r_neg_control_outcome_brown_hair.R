############################################################
## TwoSampleMR: GLP-1R instruments -> Hair colour (dark brown)
## Exposure instruments CSV: js9_2024_04_20_deng_ijs-d-24-00929r1_sdc1_READY_MR.csv
## Outcome: OpenGWAS ebi-a-GCST90029032
## Run this script line-by-line in RStudio.
############################################################
rm(list = ls())

## 0) Housekeeping: set working directory to the script location (RStudio)
if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  this_path <- rstudioapi::getActiveDocumentContext()$path
  if (nzchar(this_path)) setwd(dirname(this_path))
}
getwd()

## If you prefer an explicit path, set it here (optional):
## setwd("G:\\My Drive\\...\\some_folder")

## 1) Load packages
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
  library(readr)
})

## 2) Inputs
exposure_csv <- "js9_2024_04_20_deng_ijs-d-24-00929r1_sdc1_READY_MR.csv"
out_prefix   <- "results_prefix_glp1r_neg_control"

outcome_id   <- "ebi-a-GCST90029032"
outcome_name <- "Hair color (dark brown)"

## IMPORTANT: set OpenGWAS token before calling extract_outcome_data()
## Sys.setenv(OPENGWAS_JWT = "PASTE_YOUR_TOKEN_HERE")

## 3) Sanity check
stopifnot(file.exists(exposure_csv))

## 4) Read exposure instruments (MR-ready CSV)
exp_raw <- readr::read_csv(exposure_csv, show_col_types = FALSE)

## 5) Standardise column names if needed
## This script expects at minimum:
##   SNP, beta, se, effect_allele, other_allele, eaf, pval
## If your file already uses these exact names, this does nothing.
## If it uses common alternatives (e.g. EA/OA, BETA/SE, P), we map them below.
nm <- names(exp_raw)

rename_map <- c(
  "rsid" = "SNP",
  "RSID" = "SNP",
  "variant" = "SNP",
  "variant_id" = "SNP",
  "chr_pos" = "SNP",          # only correct if you truly stored SNP IDs here; otherwise remove
  "EA" = "effect_allele",
  "A1" = "effect_allele",
  "effectAllele" = "effect_allele",
  "OA" = "other_allele",
  "A2" = "other_allele",
  "otherAllele" = "other_allele",
  "BETA" = "beta",
  "Beta" = "beta",
  "SE" = "se",
  "Se" = "se",
  "EAF" = "eaf",
  "AF" = "eaf",
  "P" = "pval",
  "PVALUE" = "pval",
  "p_value" = "pval",
  "pvalue" = "pval",
  "N" = "samplesize",
  "n" = "samplesize"
)

# Apply renames only for columns that actually exist
to_rename <- intersect(names(rename_map), names(exp_raw))
if (length(to_rename) > 0) {
  exp_raw <- exp_raw %>% rename(!!!setNames(rename_map[to_rename], to_rename))
}

required_cols <- c("SNP","beta","se","effect_allele","other_allele","eaf","pval")
missing <- setdiff(required_cols, names(exp_raw))
if (length(missing) > 0) {
  stop(
    paste0(
      "Exposure CSV is missing required column(s): ",
      paste(missing, collapse = ", "),
      "\nPresent columns are: ", paste(names(exp_raw), collapse = ", ")
    ),
    call. = FALSE
  )
}

## 6) Coerce types (robust p-value parsing)
exp_raw <- exp_raw %>%
  mutate(
    beta = as.numeric(beta),
    se   = as.numeric(se),
    eaf  = as.numeric(eaf),
    pval = readr::parse_number(as.character(pval))
  )

bad_p <- sum(is.na(exp_raw$pval))
message("pval NAs after parsing: ", bad_p, " / ", nrow(exp_raw))

## 7) Format for TwoSampleMR (avoid samplesize_col=NULL bug)
if ("samplesize" %in% names(exp_raw)) {
  exp_dat <- TwoSampleMR::format_data(
    dat = exp_raw,
    type = "exposure",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "eaf",
    pval_col = "pval",
    samplesize_col = "samplesize"
  )
} else {
  exp_dat <- TwoSampleMR::format_data(
    dat = exp_raw,
    type = "exposure",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "eaf",
    pval_col = "pval"
  )
}

## Exposure label
exp_dat$exposure <- if ("exposure" %in% names(exp_raw)) exp_raw$exposure[1] else "GLP1R"

## Filter to valid instruments
exp_dat <- exp_dat %>%
  filter(!is.na(SNP), !is.na(beta.exposure), !is.na(se.exposure), !is.na(pval.exposure))

stopifnot(nrow(exp_dat) >= 3)
message("Exposure instruments loaded: ", nrow(exp_dat))

## 8) Pull outcome associations from OpenGWAS
out_dat <- TwoSampleMR::extract_outcome_data(
  snps = exp_dat$SNP,
  outcomes = outcome_id
)

if (is.null(out_dat) || nrow(out_dat) == 0) {
  stop("No outcome data returned for the provided SNPs and outcome ID.", call. = FALSE)
}

out_dat$outcome <- outcome_name
message("Outcome associations retrieved: ", nrow(out_dat))

## 9) Harmonise
dat <- TwoSampleMR::harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat  = out_dat,
  action = 2
)

stopifnot(nrow(dat) >= 3)
message("Instruments after harmonisation: ", nrow(dat))

## 10) MR + sensitivity
mr_res     <- TwoSampleMR::mr(dat)
het_res    <- TwoSampleMR::mr_heterogeneity(dat)
pleio_res  <- TwoSampleMR::mr_pleiotropy_test(dat)
single_res <- TwoSampleMR::mr_singlesnp(dat)
loo_res    <- TwoSampleMR::mr_leaveoneout(dat)

mr_res_or <- tryCatch(TwoSampleMR::generate_odds_ratios(mr_res), error = function(e) NULL)

## 11) Save outputs
readr::write_csv(dat,        paste0(out_prefix, "_harmonised.csv"))
readr::write_csv(mr_res,     paste0(out_prefix, "_mr_results.csv"))
readr::write_csv(het_res,    paste0(out_prefix, "_heterogeneity.csv"))
readr::write_csv(pleio_res,  paste0(out_prefix, "_pleiotropy.csv"))
readr::write_csv(single_res, paste0(out_prefix, "_singlesnp.csv"))
readr::write_csv(loo_res,    paste0(out_prefix, "_leaveoneout.csv"))
if (!is.null(mr_res_or)) readr::write_csv(mr_res_or, paste0(out_prefix, "_mr_results_OR.csv"))

## 12) Inspect results
mr_res %>% select(method, nsnp, b, se, pval) %>% print(n = Inf)

message("Done. Wrote files with prefix: ", out_prefix)
