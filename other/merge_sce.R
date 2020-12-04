library(SingleCellExperiment)
library(scran)
library(batchelor)
library(scater)
library(data.table)
library(purrr)

#####################
## Define settings ##
#####################

io <- list()

# io$sce1 <- "/Users/ricard/data/10x_gastrulation_DNMTs/processed/all_batches/SingleCellExperiment.rds"
# io$sce2 <- "/Users/ricard/data/10x_gastrulation_DNMTs/processed/sixth_batch/SingleCellExperiment.rds"
# io$metadata <-"/Users/ricard/data/10x_gastrulation_DNMTs/sample_metadata.txt.gz"
# io$outfile <- "/Users/ricard/data/10x_gastrulation_DNMTs/processed/all_batches/SingleCellExperiment.rds"

io$sce1 <- "/hps/nobackup2/research/stegle/users/ricard/10x_gastrulation_DNMTs/processed/all_batches/SingleCellExperiment.rds"
io$sce2 <- "/hps/nobackup2/research/stegle/users/ricard/10x_gastrulation_DNMTs/processed/sixth_batch/SingleCellExperiment.rds"
io$metadata <-"/hps/nobackup2/research/stegle/users/ricard/10x_gastrulation_DNMTs/sample_metadata.txt.gz"
io$outfile <- "/hps/nobackup2/research/stegle/users/ricard/10x_gastrulation_DNMTs/processed/all_batches/SingleCellExperiment.rds"

###############
## Load data ##
###############

a <- readRDS(io$sce1)
b <- readRDS(io$sce2)
sample_metadata <- fread(io$metadata)

##################
## Match genes ##
##################

all(rownames(a) == rownames(b))
# a <- a[intersect(rownames(a), rownames(b)),]
# b <- b[intersect(rownames(a), rownames(b)),]

#####################
## Combine objects ##
#####################

sce <- SingleCellExperiment(
  list(counts=Matrix::Matrix(cbind(counts(a),counts(b)),sparse=TRUE)))

# head(sizeFactors(a))

# sce <- multiBatchNorm(sce, batch=batch)
sce <- logNormCounts(sce)

#######################
## Add cell metadata ##
#######################

coldata <- sample_metadata %>%
  .[cell %in% colnames(sce)] %>% 
  setkey(cell) %>% .[colnames(sce)] %>%
  tibble::column_to_rownames("cell")
stopifnot(all(colnames(sce)==rownames(coldata)))
colData(sce) <- DataFrame(coldata)

##########
## Save ##
##########

saveRDS(sce, io$outfile)
