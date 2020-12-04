
####################
## Load libraries ##
####################

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(HDF5Array))

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
# io$outdir <- paste0(io$basedir,"/results/differential")

# Define options
# opts$test <- TRUE

#####################
## Update metadata ##
#####################

# if (opts$test) sample_metadata <- head(sample_metadata,n=100)

###############
## Load data ##
###############

# Load SingleCellExperimemt object
io$sce <- "/Users/ricard/data/10x_gastrulation_DNMTs/processed/sixth_batch/SingleCellExperiment.rds"
sce <- readRDS(io$sce)

m <- as.matrix(counts(sce))

##########################
## Create HDF5 data set ##
##########################

io$h5File <- '/Users/ricard/data/10x_gastrulation_DNMTs/processed/sixth_batch/hdf5/foo.h5'

h5createFile(io$h5File)

h5createDataset(
  file = io$h5File, 
  dataset = "counts", 
  dims = dim(m), 
  chunk = c(100,100),
  level = 0, 
  storage.mode = "integer"
)
h5write(m, file = io$h5File, name = "counts")


#################################
## Create SingleCellExperiment ##
#################################

h5.uncmp <- HDF5Array(io$h5File, name = "counts", as.sparse=TRUE)
h5read("myhdf5file.h5", "foo/S", index=list(2:3,c(1,2,4,5)))

class(h5.uncmp)

tenx.uncmp <- SingleCellExperiment(
    list(counts = h5.uncmp), 
    rowData = rowData(tenx), 
    colData = colData(tenx)
)

