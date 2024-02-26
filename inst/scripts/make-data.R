# 'muleaData' is an ExperimentHubData package for the 'mulea' comprehensive
# overrepresentation and functional enrichment analyser R package. We created
# ontologies (gene sets) in a standardised GMT (Gene Matrix Transposed) format
# for 27 different model organisms, ranging from Escherichia coli to human. The
# datasets were acquired from publicly available data sources, like gene
# expression	data from the FlyAtlas and ModEncode databases, genomic location
# data such as chromosomal bands and consecutive genes from the ensembl
# database, miRNA regulation data from the miRTarBase database, Gene ontology
# data from Gene Ontology (GO) database, pathway	data from the Pathway Commons,
# Reactome, Signalink, and Wikipathways databases, protein domain data from the
# PFAM database, transcription factor regulation	data from the ATRM, dorothEA,
# RegulonDB, TFLink, TRRUST, and Yeastract databases. The GMT files are
# available at
# https://github.com/ELTEbioinformatics/GMT_files_for_mulea/tree/main/GMT_files.

# The scripts we applied to download and process the ontologies are available at
# https://github.com/ELTEbioinformatics/GMT_files_for_mulea/tree/main/\
# scripts_to_create_GMT_files

# The script we applied to map different gene and protein identifiers, such as
# UniProt protein IDs, Entrez, Gene Symbol, and Ensembl gene IDs, is available
# at
# https://github.com/ELTEbioinformatics/GMT_files_for_mulea/tree/main/\
# scripts_to_create_GMT_files/ID_mapping_scripts

# We read and saved these ontology files using the mulea::read_gmt() and the
# saveRDS() functions using the following R script on a Linux PC:

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Saving GMTs to Rds-s ----------------------------------------------------
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# * Downloading the GMT files ----
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

# timeout should be set to larger than 60 seconds to avoid errors
options(timeout = 200)

# download a .zip file of the repository
# from the "Clone or download - Download ZIP" button
# on the GitHub repository of interest
download.file(
url =
"https://github.com/ELTEbioinformatics/GMT_files_for_mulea/archive/master.zip",
destfile = "../GMT_files_for_mulea-main.zip")

# unzip the .zip file
unzip(zipfile = "../GMT_files_for_mulea-main.zip", exdir = "..")

# delete the *.zip file
system2("rm ../GMT_files_for_mulea-main.zip")

# Set the location of the folders containing the GMT files of each taxon
GMT_dir = "../GMT_files_for_mulea-main/GMT_files"

# delete if there is a .DS_Store file in folders of GMT
system2(paste0("find ", GMT_dir, " -name '.DS_Store' -type f -delete"))

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# * Unzip zipped GMT files ----
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# GMT files bigger than 100 Mb are zipped to make it compatible with GitHub

# vector with name of the *.zip file(s) and the relative path(s) among the GMT
# files
zip_files <- system(paste0("find ", GMT_dir, " -name '*.zip'"), intern = TRUE)

# unzip the zip file(s)
lapply(zip_files, function(x) unzip(x, exdir = dirname(x)))

# delete the *.zip file(s)
system2(paste0("find ", GMT_dir, " -name '*.zip' -type f -delete"))

# delete the the __MACOSX folders
system2(paste0("find ", GMT_dir, " -name '__MACOSX' -type d -exec rm -rf {} +"))

rm(zip_files)

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# * Saving GMTs to Rds-s ----
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

destination_path <- "../rds"
# create the directory if it does not exist
if(!dir.exists(destination_path)) dir.create(destination_path)

# Get the list of files
files <- list.files(path = GMT_dir, recursive = TRUE, full.names = TRUE)

# Function to read and save each file
process_file <- function(file_path) {
    # Read GMT file
    gmt_data <- mulea::read_gmt(file_path)

    # Extract the file name without extension
    file_name <- tools::file_path_sans_ext(basename(file_path))

    # Create the destination file path for RDS
    destination_file_path <- file.path(
        destination_path,
        paste0(file_name, ".rds"))

    # Save the data in RDS format
    saveRDS(gmt_data, file = destination_file_path)

    # Print a message to indicate progress
    print("File", file_name, "processed and saved to",
            destination_file_path, "\n")
}

# Apply the function to each file
lapply(files, process_file)

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Session info
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

sessionInfo()
# R version 4.3.2 (2023-10-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 22.04.3 LTS
#
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0
#
# locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=hu_HU.UTF-8
# LC_COLLATE=en_US.UTF-8
# [5] LC_MONETARY=hu_HU.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=hu_HU.UTF-8
# LC_NAME=C
# [9] LC_ADDRESS=C               LC_TELEPHONE=C
# LC_MEASUREMENT=hu_HU.UTF-8 LC_IDENTIFICATION=C
#
# time zone: Europe/Budapest
# tzcode source: system (glibc)
#
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base
#
# loaded via a namespace (and not attached):
# [1] compiler_4.3.2 tools_4.3.2
