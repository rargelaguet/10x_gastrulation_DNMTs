
####################
## Load libraries ##
####################

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/10x_gastrulation_DNMTs/settings.R")
} else {
  stop("Computer not recognised")
}

# Define I/O
io$outdir <- paste0(io$basedir,"/results/differential")

# Define options
opts$test <- TRUE

#####################
## Update metadata ##
#####################

if (opts$test) sample_metadata <- head(sample_metadata,n=100)

###############
## Load data ##
###############

# Load SingleCellExperimemt object
sce <- readRDS(io$sce)
