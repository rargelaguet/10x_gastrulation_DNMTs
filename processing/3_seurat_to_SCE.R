library(scater)
library(scran)

#####################
## Define settings ##
#####################

# Load default settings
# source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
source("/homes/ricard/10x_gastrulation_DNMTs/settings.R")

# Define I/O
# io$srat <- "/Users/ricard/data/10x_gastrulation_DNMTs/processed/first_batch/srat.rds"
io$outfile <- io$sce

# io$srat <- "/hps/nobackup2/research/stegle/users/ricard/10x_gastrulation_DNMTs/processed/first_batch/srat.rds"
# io$outfile <- "/hps/nobackup2/research/stegle/users/ricard/10x_gastrulation_DNMTs/processed/first_batch/SingleCellExperiment.rds"
# io$metadata <- "/hps/nobackup2/research/stegle/users/ricard/10x_gastrulation_DNMTs/processed/first_batch/sample_metadata.txt.gz"

###############
## Load data ##
###############

# Load sample metadata
sample_metadata <- fread(io$metadata)

# Load seurat
srat <- readRDS(io$srat)

#####################################
## Convert to SingleCellExperiment ##
#####################################

sce <- as.SingleCellExperiment(srat)

# Add metadata
# stopifnot(sample_metadata$cell%in%colnames(sce))
# stopifnot(colnames(sce)%in%sample_metadata$cell)
sample_metadata <- sample_metadata %>% .[cell%in%colnames(sce)] %>% setkey(cell) %>% .[colnames(sce)]
stopifnot(sample_metadata$cell == colnames(sce))
colData(sce) <- sample_metadata %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
  .[colnames(sce),] %>% DataFrame()

# Normalise using scran
# clusts = as.numeric(quickCluster(sce, method = "igraph", min.size = 100, BPPARAM = mcparam))
clusts = as.numeric(quickCluster(sce))
min.clust = min(table(clusts))/2
new_sizes = c(floor(min.clust/3), floor(min.clust/2), floor(min.clust))
# sce = computeSumFactors(sce, clusters = clusts, sizes = new_sizes, max.cluster.size = 3000)
sce = computeSumFactors(sce)

to.plot <- data.frame(X = Matrix::colSums(counts(sce)), Y = sizeFactors(sce))
ggplot(to.plot, mapping = aes(x = X, y = Y)) +
  geom_point() +
  labs(x = "Number of UMIs", y = "Size Factor") +
  theme_classic()

sce <- logNormCounts(sce)

saveRDS(sce, io$outfile)
