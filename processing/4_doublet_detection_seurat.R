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
p$add_argument('--samples',                 type="character",   nargs='+',     help='Sample(s)')
p$add_argument('--doublet_score_threshold',  type="double",   help='Doublet score threshold')
p$add_argument('--test',                    action = "store_true",             help='Testing mode')
p$add_argument('--outfile',                  type="character",                  help='Output file')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################


## START TEST ##
args <- list()
args$seurat <- file.path(io$basedir,"processed_new/Seurat.rds")
args$metadata <- file.path(io$basedir,"results_new/qc/sample_metadata_after_qc.txt.gz")
args$samples <- opts$samples[1]
args$hybrid_score_threshold <- 1.50
args$test <- FALSE
## END TEST ##

if (isTRUE(args$test)) print("Test mode activated...")

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_rnaQC==TRUE & sample%in%args$samples]
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

#############################
## Calculate doublet score ##
#############################

seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat)

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list <- paramSweep_v3(seurat, PCs = 1:10, sct = FALSE)
sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
bcmvn_kidney <- find.pK(sweep.stats_kidney)


# Call doublets
dt[,doublet_call:=hybrid_score>args$hybrid_score_threshold]
table(dt$doublet_call)

##########
## Save ##
##########

fwrite(dt, args$outfile, sep="\t", na="NA", quote=F)


