suppressPackageStartupMessages(library(DoubletFinder))
suppressPackageStartupMessages(library(Seurat))

here::i_am("processing/4_doublet_detection.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--seurat',         type="character",                            help='SingleCellExperiment file')
p$add_argument('--metadata',    type="character",                            help='metadata file')
p$add_argument('--sample',                 type="character",     help='Sample(s)')
p$add_argument('--number_doublets',  type="integer",   help='Number of doublets')
p$add_argument('--outfile',                  type="character",                  help='Output file')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
# args <- list()
# args$seurat <- file.path(io$basedir,"processed_new/seurat.rds")
# args$metadata <- file.path(io$basedir,"results_new/qc/sample_metadata_after_qc.txt.gz")
# args$sample <- opts$samples[1]
# args$number_doublets <- 100
## END TEST ##

if (isTRUE(args$test)) print("Test mode activated...")

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_rnaQC==TRUE & sample==args$sample]
table(sample_metadata$sample)

###############
## Load data ##
###############

# Load RNA expression data as SingleCellExperiment object
seurat <- load_Seurat(
  file = args$seurat, 
  cells = sample_metadata$cell,
  normalise = TRUE, scale = FALSE, 
  remove_non_expressed_genes = TRUE
)
dim(seurat)

# Update sample metadata
foo <- sample_metadata %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
  .[colnames(seurat),]
stopifnot(colnames(seurat) == rownames(foo))
seurat@meta.data <- foo

##############################
## Dimensionality reduction ##
##############################

seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat, npcs=25)

#############################
## Calculate doublet score ##
#############################

# Performs pN-pK parameter sweeps on a 10,000-cell subset
# Parameter ranges tested: pN = 0.05-0.3, pK = 0.0005-0.3.
# Outputs a list of pANN vectors for every pN and pK combination. Output also contains pANN information for artificial doublets.
sweep.res.list <- paramSweep_v3(seurat, PCs = 1:25, sct = FALSE)

# Summarizes results from doubletFinder_ParamSweep, computing the bimodality coefficient across pN and pK parameter space
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)

# Computes and visualizes the mean-variance normalized bimodality coefficient (BCmvn) score for each pK value
# Optimal pK can be manually discerned as maxima in BCmvn distributions
# Returns a Dataframe of mean BC, BC variance, and BCmvn scores for each pK value.
tmp <- find.pK(sweep.stats)

# Fetch optimal pK
pK.optim <- tmp[which.max(tmp$BCmetric),"pK"] %>% as.character %>% as.numeric

# Homotypic Doublet Proportion Estimate 
# Leverages user-provided cell annotations to model the proportion of homotypic doublets
# homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
# nExp_poi <- round(0.075*nrow(seu_kidney@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
# nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder with optimal hyperparameters
# seurat <- doubletFinder_v3(seurat, 1:25, pN = 0.25, pK = pK.optim, nExp = as.integer(0.01*ncol(seurat)), reuse.pANN = FALSE)
seurat <- doubletFinder_v3(seurat, 1:25, pN = 0.25, pK = pK.optim, nExp = args$number_doublets, reuse.pANN = FALSE)

# Call doublets
dt <- data.table(
  cell = colnames(seurat),
  # sample = args$sample,
  doublet_score = seurat@meta.data[[grep("pANN",colnames(seurat@meta.data),value=T)]] %>% round(3),
  doublet_call = seurat@meta.data[[grep("classifications",colnames(seurat@meta.data),value=T)]]=="Doublet"
)

print("Number of doublets:")
print(table(dt$doublet_call))

##########
## Save ##
##########

fwrite(dt, args$outfile, sep="\t", na="NA", quote=F)

