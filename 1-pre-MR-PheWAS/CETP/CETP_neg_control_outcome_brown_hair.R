############################################################
## TwoSampleMR: CETP instruments -> Hair colour (brown)
## Outcome: OpenGWAS ebi-a-GCST90029032
## Run this script line-by-line in RStudio.
## Assumes your working directory contains:
##   - this script
##   - CETP_blauw_2018_DOI_10.1161_CIRCGEN.117.00203.csv
############################################################
rm(list = ls())

## 0) Housekeeping: set working directory to the script location (RStudio)
if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  this_path <- rstudioapi::getActiveDocumentContext()$path
  if (nzchar(this_path)) setwd(dirname(this_path))
}
getwd()

## If you prefer an explicit path, keep this (otherwise remove it):
setwd("G:\\My Drive\\Documents\\0Oxford_main\\ARD paper\\4_paper_stuff\\ad hoc analysis\\CETP")

## 1) Load packages
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
  library(readr)
})

## 2) Inputs
exposure_csv <- "CETP_blauw_2018_DOI_10.1161_CIRCGEN.117.00203.csv"

# Use a distinct prefix so you do not overwrite earlier outputs
out_prefix <- "CETP_neg_control_hair_brown"

outcome_id   <- "ebi-a-GCST90029032"
outcome_name <- "Hair colour: brown"

## IMPORTANT: set OpenGWAS token before calling extract_outcome_data()
## Either set it here (session-only)...
## Sys.setenv(OPENGWAS_JWT = "PASTE_YOUR_TOKEN_HERE")

## ...or ensure it is set in .Renviron and restart R.

## 3) Sanity check
stopifnot(file.exists(exposure_csv))

## 4) Read exposure instruments (pre-QC'd)
exp_raw <- readr::read_csv(exposure_csv, show_col_types = FALSE)

required_cols <- c("SNP","beta","se","effect_allele","other_allele","eaf","pval")
missing <- setdiff(required_cols, names(exp_raw))
if (length(missing) > 0) {
  stop(
    paste0("Exposure CSV is missing required column(s): ", paste(missing, collapse = ", ")),
    call. = FALSE
  )
}

## 5) Format for TwoSampleMR (robust to missing samplesize + pval parsing)
exp_raw <- exp_raw %>%
  mutate(
    beta = as.numeric(beta),
    se   = as.numeric(se),
    eaf  = as.numeric(eaf),
    pval = readr::parse_number(as.character(pval))
  )

bad_p <- sum(is.na(exp_raw$pval))
message("pval NAs after parsing: ", bad_p, " / ", nrow(exp_raw))

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

## Optional exposure label
exp_dat$exposure <- if ("exposure" %in% names(exp_raw)) exp_raw$exposure[1] else "CETP"

## Keep only valid rows (and drop missing p-values too)
exp_dat <- exp_dat %>%
  filter(!is.na(SNP), !is.na(beta.exposure), !is.na(se.exposure), !is.na(pval.exposure))

stopifnot(nrow(exp_dat) >= 3)
message("Exposure instruments loaded: ", nrow(exp_dat))

## 6) Pull outcome associations from OpenGWAS
out_dat <- TwoSampleMR::extract_outcome_data(
  snps = exp_dat$SNP,
  outcomes = outcome_id
)

if (is.null(out_dat) || nrow(out_dat) == 0) {
  stop("No outcome data returned for the provided SNPs and outcome ID.", call. = FALSE)
}

out_dat$outcome <- outcome_name
message("Outcome associations retrieved: ", nrow(out_dat))

## 7) Harmonise
dat <- TwoSampleMR::harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat  = out_dat,
  action = 2
)

stopifnot(nrow(dat) >= 3)
message("Instruments after harmonisation: ", nrow(dat))

## 8) MR + sensitivity
mr_res     <- TwoSampleMR::mr(dat)
het_res    <- TwoSampleMR::mr_heterogeneity(dat)
pleio_res  <- TwoSampleMR::mr_pleiotropy_test(dat)
single_res <- TwoSampleMR::mr_singlesnp(dat)
loo_res    <- TwoSampleMR::mr_leaveoneout(dat)

## Hair colour is typically binary/ordinal-coded in GWAS; OR conversion may or may not be meaningful.
## Keep as optional.
mr_res_or <- tryCatch(TwoSampleMR::generate_odds_ratios(mr_res), error = function(e) NULL)

## 9) Save outputs
readr::write_csv(dat,        paste0(out_prefix, "_harmonised.csv"))
readr::write_csv(mr_res,     paste0(out_prefix, "_mr_results.csv"))
readr::write_csv(het_res,    paste0(out_prefix, "_heterogeneity.csv"))
readr::write_csv(pleio_res,  paste0(out_prefix, "_pleiotropy.csv"))
readr::write_csv(single_res, paste0(out_prefix, "_singlesnp.csv"))
readr::write_csv(loo_res,    paste0(out_prefix, "_leaveoneout.csv"))
if (!is.null(mr_res_or)) readr::write_csv(mr_res_or, paste0(out_prefix, "_mr_results_OR.csv"))

## 10) Inspect results
mr_res %>% select(method, nsnp, b, se, pval) %>% print(n = Inf)

message("Done. Wrote files with prefix: ", out_prefix)
