suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(scran))

#####################
## Define settings ##
#####################

# Load default settings
source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
# source("/homes/ricard/10x_gastrulation_DNMTs/settings.R")

# Define I/O
io$metadata <- "/Users/ricard/data/10x_gastrulation_DNMTs/processed/sixth_batch/metadata.txt.gz"
io$seurat <- "/Users/ricard/data/10x_gastrulation_DNMTs/processed/sixth_batch/seurat.rds"
io$sce <- "/Users/ricard/data/10x_gastrulation_DNMTs/processed/sixth_batch/SingleCellExperiment.rds"
io$outfile <- io$sce

# Define options
opts$test <- FALSE

###############
## Load data ##
###############

# Load sample metadata
sample_metadata <- fread(io$metadata) %>% .[pass_QC==TRUE]

if (opts$test) {
	sample_metadata <- sample_metadata %>% head(n=1000)	
}
# Load seurat
seurat <- readRDS(io$seurat)[,sample_metadata$cell]

#####################################
## Convert to SingleCellExperiment ##
#####################################

sce <- as.SingleCellExperiment(seurat)

# Add metadata
# stopifnot(sample_metadata$cell%in%colnames(sce))
# stopifnot(colnames(sce)%in%sample_metadata$cell)
sample_metadata <- sample_metadata %>% .[cell%in%colnames(sce)] %>% setkey(cell) %>% .[colnames(sce)]
stopifnot(sample_metadata$cell == colnames(sce))
colData(sce) <- sample_metadata %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
  .[colnames(sce),] %>% DataFrame()

##########################
## Compute size factors ##
##########################

# clusts = as.numeric(quickCluster(sce, method = "igraph", min.size = 100, BPPARAM = mcparam))
# clusts = as.numeric(quickCluster(sce))
# min.clust = min(table(clusts))/2
# new_sizes = c(floor(min.clust/3), floor(min.clust/2), floor(min.clust))
# sce = computeSumFactors(sce, clusters = clusts, sizes = new_sizes, max.cluster.size = 3000)
sce = computeSumFactors(sce)

###########################
## Normalise using scran ##
###########################

sce <- logNormCounts(sce)

##########
## Plot ##
##########

# to.plot <- data.frame(X = Matrix::colSums(counts(sce)), Y = sizeFactors(sce))
# ggplot(to.plot, mapping = aes(x = X, y = Y)) +
#   geom_point() +
#   labs(x = "Number of UMIs", y = "Size Factor") +
#   theme_classic()

##########
## Save ##
##########

saveRDS(sce, io$outfile)
