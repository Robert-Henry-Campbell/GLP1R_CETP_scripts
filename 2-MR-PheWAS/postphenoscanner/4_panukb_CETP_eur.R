#panukb run

rm(list = ls())

setwd('G:\\My Drive\\Documents\\0Oxford_main\\ARD paper\\ardmr')
#pull
#system("git fetch origin")
#system("git pull origin main")


#save chagen and make them avalible in local session
devtools::load_all()
Sys.setenv(ARDMR_CACHE_DIR = 'C:\\Users\\Robert\\Downloads\\ardmr')



# snps
cetp_snps <- read.csv("G:\\My Drive\\Documents\\0Oxford_main\\ARD paper\\3_ARD_MR\\exposures to run on\\CETP\\CETP_blauw 2018 DOI_10.1161 CIRCGEN.117.00203.csv")


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

