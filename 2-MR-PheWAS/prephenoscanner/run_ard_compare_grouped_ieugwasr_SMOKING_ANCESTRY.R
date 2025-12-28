#run_ard_compare_grouped_ieugwasr

rm(list = ls())

script_dir <- dirname(this.path::this.path())
setwd(normalizePath(file.path(script_dir, "..", "..")))
#set cache dir
Sys.setenv(ARDMR_CACHE_DIR = '3-output')
library(ardmr)

#set token
Sys.setenv(OPENGWAS_JWT = "YOUR TOKEN GOES HERE")



#set csv path
csv_path = '2-MR-PheWAS/exposures/ieugwasr pre phenoscanner/ieugwasr_exposures_23.9.25_smoking_only.csv'

#call the ard_compare_grouped_ieugwasr
result <- run_ieugwasr_ard_compare(
  csv_path = csv_path,
  cache_dir = ardmr_cache_dir(),   # or your custom path
  prompt_for_units = TRUE          # will ask you for any missing units
)
