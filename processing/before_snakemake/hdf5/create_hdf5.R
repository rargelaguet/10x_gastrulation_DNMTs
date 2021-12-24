
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
h5write(m, io$h5File, name = "counts")


# print(object.size(m), units="auto")

#################
## HDF5 matrix ##
#################

# Does not read it into memory
m.hdf5 <- HDF5Array(io$h5File, name = "counts", as.sparse=TRUE)
class(m.hdf5)
dim(m.hdf5)
print(object.size(m.hdf5), units="auto")
is(m.hdf5, "DelayedArray")

# Reads it into memory
m.hdf5 <- h5read(io$h5File, "counts")
class(m.hdf5)
dim(m.hdf5)
print(object.size(m.hdf5), units="auto")

colnames(m.hdf5) <- colnames(sce)
rownames(m.hdf5) <- rownames(sce)

################################
## Create DelayedArray matrix ##
################################

library(DelayedArray)

foo <- m.hdf5[,1:100]

#################################
## Create SingleCellExperiment ##
#################################


sce.hdf5 <- SingleCellExperiment(
    list(counts = m.hdf5), 
    rowData = rowData(sce), 
    colData = colData(sce)
)
print(object.size(sce.hdf5), units="auto")

