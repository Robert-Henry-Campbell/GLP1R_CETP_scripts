#panukb run

rm(list = ls())

script_dir <- dirname(this.path::this.path())
setwd(normalizePath(file.path(script_dir, "..", "..")))
Sys.setenv(ARDMR_CACHE_DIR = '3-output')

library(ardmr)

Sys.setenv(OPENGWAS_JWT = "YOUR TOKEN GOES HERE")

# snps
cetp_snps <- read.csv("2-MR-PheWAS/exposures/CETP/CETP_blauw 2018 DOI_10.1161 CIRCGEN.117.00203.csv")


CETP_run <- run_phenome_mr(exposure = 'CETP',
                           exposure_units = 'µg/mL CETP serum level',
                           exposure_snps = cetp_snps,
                           ancestry = 'EUR',
                           sex = 'both',
                           Multiple_testing_correction = 'BH',
                           confirm = 'yes',
                           force_refresh = TRUE,
                           sensitivity_pass_min = 5
)
