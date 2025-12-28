rm(list = ls())

script_dir <- dirname(this.path::this.path())
setwd(normalizePath(file.path(script_dir, "..", "..")))
Sys.setenv(ARDMR_CACHE_DIR = '3-output')

library(ardmr)

#set token
Sys.setenv(OPENGWAS_JWT = "YOUR TOKEN GOES HERE")

exposure_snps <- read.csv("2-MR-PheWAS/exposures/glp1/js9_2024_04_20_deng_ijs-d-24-00929r1_sdc1_READY_MR.csv")


GLP1_run <- run_phenome_mr(exposure = 'GLP1R',
                       exposure_units = 'SD GLP1 expression',
                       exposure_snps = exposure_snps,
                       ancestry = 'EUR',
                       sex = 'both',
                       Multiple_testing_correction = 'BH',
                       confirm = 'yes',
                       force_refresh = TRUE,
                       sensitivity_pass_min = 6
)

