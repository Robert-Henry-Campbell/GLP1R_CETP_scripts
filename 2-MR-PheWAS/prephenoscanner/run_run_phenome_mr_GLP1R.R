#pull
system("git fetch origin")
system("git pull origin main")

setwd('G:\\My Drive\\Documents\\0Oxford_main\\ARD paper\\ardmr')
#save chagen and make them avalible in local session
devtools::load_all()
rm(list = ls())
#example snps
exposure_snps <- read.csv("G:\\My Drive\\Documents\\0Oxford_main\\ARD paper\\3_ARD_MR\\exposures to run on\\glp1\\js9_2024_04_20_deng_ijs-d-24-00929r1_sdc1_READY_MR.csv")
Sys.setenv(ARDMR_CACHE_DIR = 'C:\\Users\\Robert\\Downloads\\ardmr')

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


