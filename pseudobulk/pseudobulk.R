# suppressPackageStartupMessages(library(Seurat))

source(here::here("settings.R"))
source(here::here("utils.R"))
source(here::here("pseudobulk/utils.R"))

here::i_am("pseudobulk/pseudobulk.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sce',         type="character",    help='SingleCellExperiment object')
p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--group_by',    type="character",    help='Metadata column to group cells by')
p$add_argument('--normalisation_method',    type="character",    help='Metadata column to group cells by')
p$add_argument('--outdir',      type="character",    help='Output directory')

args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
# args$metadata <- file.path(io$basedir,"results/mapping/sample_metadata_after_mapping.txt.gz")
# args$sce <- file.path(io$basedir,"processed_all/SingleCellExperiment.rds")
# args$group_by <- "class_sample_celltype"
# args$normalisation_method <- "cpm"
# args$outdir <- file.path(io$basedir,"results/pseudobulk")
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings=F)

# Options
# opts$rename_celltypes <- c(
#   "Erythroid3" = "Erythroid",
#   "Erythroid2" = "Erythroid",
#   "Erythroid1" = "Erythroid",
#   "Blood_progenitors_1" = "Blood_progenitors",
#   "Blood_progenitors_2" = "Blood_progenitors",
#   "Anterior_Primitive_Streak" = "Primitive_Streak"
# )


###############
## Load data ##
###############

# Load cell metadata
sample_metadata <- fread(args$metadata) %>%
  .[,dataset:=ifelse(grepl("Grosswendt",sample),"CRISPR","KO")] %>%
  # .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,opts$rename_celltypes)] %>%
  .[,class_celltype:=sprintf("%s-%s",class,celltype.mapped)] %>%
  .[,class_celltype_dataset:=sprintf("%s-%s-%s",class,celltype.mapped,dataset)] %>%
  .[,class_sample_celltype:=sprintf("%s-%s-%s",class,sample,celltype.mapped)] %>%
  # .[,class_celltype_dataset:=sprintf("%s-%s-%s",class,celltype.mapped,dataset)] %>%
  .[pass_rnaQC==TRUE & !is.na(eval(as.name(args$group_by)))]

table(sample_metadata[[args$group_by]])

# Load SingleCellExperiment
sce <- load_SingleCellExperiment(args$sce, cells=sample_metadata$cell)
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame


################
## Pseudobulk ##
################

sce_pseudobulk <- pseudobulk_sce_fn(
  x = sce,
  assay = "counts",
  by = args$group_by,
  fun = "sum",
  scale = FALSE # Should pseudo-bulks be scaled with the effective library size & multiplied by 1M?
)

assayNames(sce_pseudobulk) <- "counts"

###################
## Normalisation ##
###################

if (args$normalisation_method=="deseq2") {
  
  suppressPackageStartupMessages(library(DESeq2))
  dds <- DESeqDataSet(sce_pseudobulk, design=~1)
  dds <- varianceStabilizingTransformation(dds)
  logcounts(sce_pseudobulk) <- assay(dds)
  
} else if (args$normalisation_method=="cpm") {
  
  logcounts(sce_pseudobulk) <- log2(1e6*(sweep(counts(sce_pseudobulk),2,colSums(counts(sce_pseudobulk)),"/"))+1)
  # logcounts(sce_pseudobulk) <- edgeR::cpm(counts(sce_pseudobulk), log=TRUE, prior.count = 1)
  
} else {
  stop("Normalisation method not recognised")
}

# Save
saveRDS(sce_pseudobulk, file.path(args$outdir,sprintf("SingleCellExperiment_pseudobulk_%s.rds",args$group_by)))

########################################################
## Create Seurat object from the SingleCellExperiment ##
########################################################

# sce_pseudobulk@colData$sample <- rownames(sce_pseudobulk@colData)  # at least one metadata column is needed
# seurat_pseudobulk <- as.Seurat(sce_pseudobulk)
# seurat_pseudobulk <- RenameAssays(seurat_pseudobulk, originalexp="RNA")

# Save
# saveRDS(seurat_pseudobulk, file.path(args$outdir,sprintf("Seurat_pseudobulk_%s.rds",args$group_by)))

#####################
## Save statistics ##
#####################

stats.dt <- data.table(table(sce[[args$group_by]])) %>% setnames(c("sample","N"))
fwrite(stats.dt, file.path(args$outdir,sprintf("pseudobulk_stats_%s.txt.gz",args$group_by)), sep="\t")
