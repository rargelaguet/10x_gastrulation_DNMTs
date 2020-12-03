suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(argparse))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--batches',   type="character",   nargs='+',  help='Query batch(es)')
p$add_argument('--test',      action = "store_true",          help='Testing mode')
p$add_argument('--features',  type="integer",  default=1000,                help='Number of cores')
p$add_argument('--npcs',  type="integer", default=30,                 help='Number of cores')
p$add_argument('--outdir',    type="character",               help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$batches <- c(
#   "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001"
# )
# 
# args$features <- 1000
# args$npcs <- 30
# args$test <- TRUE
# args$outdir <- paste0(io$basedir,"/results/dimensionality_reduction")
## END TEST ##

################
## Define I/O ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
} else {
  source("/homes/ricard/10x_gastrulation_DNMTs/settings.R")
}

if (isTRUE(args$test)) print("Test mode activated...")

############################
## Update sample metadata ##
############################

sample_metadata <- fread(io$metadata) %>%
  .[pass_QC==TRUE & batch%in%args$batches]
table(sample_metadata$batch)

###############
## Load data ##
###############

# Load RNA expression data as SingleCellExperiment object
sce <- readRDS(io$sce)[,sample_metadata$cell]

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

#######################
## Feature selection ##
#######################

decomp <- modelGeneVar(sce)
decomp <- decomp[decomp$mean > 0.01,]

hvgs <- decomp[order(decomp$FDR),] %>% head(n=args$features) %>% rownames
# hvgs <- rownames(decomp)[decomp$p.value <= 0.05]

# Subset SingleCellExperiment
sce_filt <- sce[hvgs,]
dim(sce_filt)

##############################
## Dimensionality reduction ##
##############################

data <- scale(t(logcounts(sce_filt)), center = T, scale = F)

# PCA
npcs <- args$npcs
reducedDim(sce_filt, "PCA") <- irlba::prcomp_irlba(data, n=npcs)$x#[,1:npcs]

# Filter PCA solution
# reducedDim(sce_filt, "PCA") <- reducedDim(sce_filt, "PCA")[,-3]

# UMAP
set.seed(42)
sce_filt <- runUMAP(sce_filt, dimred="PCA", n_neighbors = 30, min_dist = 0.3)

##########
## Plot ##
##########

to.plot <- reducedDim(sce_filt,"UMAP") %>% as.data.table %>% 
  .[,cell:=colnames(sce_filt)] %>%
  merge(sample_metadata, by="cell")

p <- ggplot(to.plot, aes(x=V1, y=V2, fill=celltype.mapped)) +
  geom_point(size=2, shape=21, stroke=0.05) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position="none",
    legend.title=element_blank()
  )

pdf(sprintf("%s/%s_umap_%d_%d.pdf",args$outdir, paste(args$batches,collapse="-"), args$features,args$npcs), width=7, height=5)
print(p)
dev.off()