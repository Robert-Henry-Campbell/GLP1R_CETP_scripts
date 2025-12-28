To reproduce the environment:

BASH:
###
git clone https://github.com/Robert-Henry-Campbell/GLP1R_CETP_scripts.git
cd GLP1R_CETP_scripts
R
###

R:
###
install.packages("renv")   # if needed
renv::restore()            # reads renv.lock and installs the same package versions
###

after this, each script should run normally. Results for the pre-MR-PheWAS will be written to their local directory. Results from MR-PheWAS will be written to /3-output/. 

Notes: 

To run the positive and negative outcome control MR, you will need:
1. an internet connection
2. an openGWAS API token (https://api.opengwas.io/api/)


To run the both-sex MR-PheWAS analyses, you will need:
1. an internet connection
2. an openGWAS API token (https://api.opengwas.io/api/)
3. approximately 10gb of available space for downloaded files and results

To run the sex-specific MR-PheWAS analyses, you will need:
1. a high speed internet connection
2. an openGWAS API token (https://api.opengwas.io/api/)
3. approximately 800gb of available space for downloaded files and results. Please allow several days for the download of the NEALE lab GWAS summary statistics. 

