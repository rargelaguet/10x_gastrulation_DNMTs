suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(argparse))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sce',         type="character",                            help='SingleCellExperiment file')
p$add_argument('--metadata',    type="character",                            help='metadata file')
p$add_argument('--samples',     type="character",                nargs='+',  help='Query batch(es)')
p$add_argument('--features',    type="integer",    default=1000,             help='Number of cores')
p$add_argument('--npcs',        type="integer",    default=30,               help='Number of cores')
p$add_argument('--n_neighbors', type="integer",    default=30,               help='(UMAP) Number of neighbours')
p$add_argument('--min_dist',    type="double",     default=0.3,              help='(UMAP) Minimum distance')
p$add_argument('--colour_by',   type="character",  default="celltype.mapped",  nargs='+',help='(UMAP) Minimum distance')
p$add_argument('--outdir',      type="character",                                  help='Output directory')
# p$add_argument('--test',      action = "store_true",                       help='Testing mode')

args <- p$parse_args(commandArgs(TRUE))

################
## Define I/O ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
} else {
  source("/homes/ricard/10x_gastrulation_DNMTs/settings.R")
}


## START TEST ##
args$samples <- c("SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001", "SIGAB6_E85_3_Dnmt3aWT_Dnmt3b_WT_L002")
args$sce <- "/Users/ricard/data/10x_gastrulation_DNMTs/processed/all_batches/test/SingleCellExperiment.rds" # io$sce
args$metadata <- io$metadata
args$features <- 1000
args$npcs <- 30
args$test <- TRUE
args$colour_by <- c("celltype.mapped","batch")
args$outdir <- paste0(io$basedir,"/results/dimensionality_reduction")
## END TEST ##

# if (isTRUE(args$test)) print("Test mode activated...")

############################
## Update sample metadata ##
############################

sample_metadata <- fread(args$metadata) %>%
  .[pass_QC==TRUE & batch%in%args$samples]
table(sample_metadata$batch)

stopifnot(args$colour_by %in% colnames(sample_metadata))

###############
## Load data ##
###############

# Load RNA expression data as SingleCellExperiment object
# sce <- readRDS(io$sce)[,sample_metadata$cell]
sce <- load_SingleCellExperiment(args$sce, cells=sample_metadata$cell)

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

# PCA
# data <- scale(t(logcounts(sce_filt)), center = T, scale = F)
# reducedDim(sce_filt, "PCA") <- irlba::prcomp_irlba(data, n=args$npcs)$x#[,1:npcs]
sce_filt <- runPCA(sce_filt, ncomponents = args$npcs, ntop=args$features)

# Filter PCA solution
# reducedDim(sce_filt, "PCA") <- reducedDim(sce_filt, "PCA")[,-3]

# UMAP
set.seed(42)
sce_filt <- runUMAP(sce_filt, dimred="PCA", n_neighbors = args$n_neighbors, min_dist = args$min_dist)

##########
## Plot ##
##########

to.plot <- reducedDim(sce_filt,"UMAP") %>% as.data.table %>% 
  .[,cell:=colnames(sce_filt)] %>%
  merge(sample_metadata, by="cell")

for (i in args$colour_by) {
  
  # Plot
  p <- ggplot(to.plot, aes_string(x="V1", y="V2", fill=i)) +
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
  
  # Save
  outfile <- sprintf("%s/%s_umap_%d_%d_%d_%s_%s.pdf",args$outdir, paste(args$samples,collapse="-"), args$features,args$npcs,args$n_neighbors,args$min_dist,i)
  pdf(outfile, width=7, height=5)
  print(p)
  dev.off()
}

