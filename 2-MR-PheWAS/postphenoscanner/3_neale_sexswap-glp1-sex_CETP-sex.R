#run_ard_compare_grouped_ieugwasr
rm(list = ls())

script_dir <- dirname(this.path::this.path())
setwd(normalizePath(file.path(script_dir, "..", "..")))
Sys.setenv(ARDMR_CACHE_DIR = '3-output')
#set cache dir
Sys.setenv(ARDMR_CACHE_DIR = '3-output')

library(ardmr)

#set token
Sys.setenv(OPENGWAS_JWT = "YOUR TOKEN GOES HERE")

#set csv paths and snps
glp1r_snps <-  read.csv("2-MR-PheWAS/exposures/glp1/js9_2024_04_20_deng_ijs-d-24-00929r1_sdc1_READY_MR.csv")
cetp_male_snps <- read.csv("2-MR-PheWAS/exposures/CETP/CETP_male_MR_ready.csv")
cetp_female_snps <- read.csv("2-MR-PheWAS/exposures/CETP/CETP_female_MR_ready.csv")
cetp_both_snps <- read.csv("2-MR-PheWAS/exposures/CETP/CETP_blauw 2018 DOI_10.1161 CIRCGEN.117.00203.csv")


####GLP1 by SEX
glp1_groups <- list(
  list(sex = "male",   ancestry = "EUR", exposure_snps = glp1r_snps),
  list(sex = "female",   ancestry = "EUR", exposure_snps = glp1r_snps)
)
glp1_sex <- ard_compare(
  exposure = "GLP1R",                          # or your exposure label
  exposure_units = 'SD GLP1 expression',
  groups = glp1_groups,
  Multiple_testing_correction = "BH",                  # or "bonferroni"
  sensitivity_enabled = c(
    "egger_intercept","egger_slope_agreement",
    "weighted_median","weighted_mode",
    "steiger_direction","leave_one_out",
    "ivw_Q","ivw_I2"
  ),
  sensitivity_pass_min = 6,
  scatterplot = TRUE,
  snpforestplot = TRUE,
  leaveoneoutplot = TRUE,
  cache_dir = ardmr_cache_dir(),                       # or a custom path
  logfile = NULL,
  verbose = TRUE,
  confirm = "no",                                      # skip interactive prompts
  force_refresh = TRUE
)

####CETP by SEX
cetp_groups <- list(
  list(sex = "male",   ancestry = "EUR", exposure_snps = cetp_male_snps),
  list(sex = "female",   ancestry = "EUR", exposure_snps = cetp_female_snps)
)
CETP_sex <- ard_compare(
  exposure = "CETP",                          # or your exposure label
  exposure_units = 'SD CETP expression',
  groups = cetp_groups,
  Multiple_testing_correction = "BH",                  # or "bonferroni"
  sensitivity_enabled = c(
    "egger_intercept","egger_slope_agreement",
    "weighted_median","weighted_mode",
    "steiger_direction","leave_one_out",
    "ivw_Q","ivw_I2"
  ),
  sensitivity_pass_min = 5,
  scatterplot = TRUE,
  snpforestplot = TRUE,
  leaveoneoutplot = TRUE,
  cache_dir = ardmr_cache_dir(),                       # or a custom path
  logfile = NULL,
  verbose = TRUE,
  confirm = "no",                                      # skip interactive prompts
  force_refresh = TRUE
)
