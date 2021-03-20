suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(argparse))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sce',             type="character",                               help='SingleCellExperiment file')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--samples',         type="character",                nargs='+',     help='Samples')
p$add_argument('--features',        type="integer",    default=1000,                help='Number of features')
p$add_argument('--npcs',            type="integer",    default=30,                  help='Number of PCs')
p$add_argument('--vars.to.regress', type="character",                nargs='+',     help='Metadata columns to regress out')
p$add_argument('--batch.correction',type="character",                               help='Metadata column to apply batch correction on')
p$add_argument('--n_neighbors',     type="integer",    default=30,   nargs='+',     help='(UMAP) Number of neighbours')
p$add_argument('--min_dist',        type="double",     default=0.3,  nargs='+',     help='(UMAP) Minimum distance')
p$add_argument('--colour_by',       type="character",  default="celltype.mapped",  nargs='+',  help='Metadata columns to colour the UMAP by')
p$add_argument('--seed',            type="integer",    default=42,                  help='Random seed')
p$add_argument('--outdir',          type="character",                               help='Output directory')
# p$add_argument('--test',      action = "store_true",                       help='Testing mode')

args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
  source("/Users/ricard/10x_gastrulation_DNMTs/utils.R")
} else {
  source("/homes/ricard/10x_gastrulation_DNMTs/settings.R")
  source("/homes/ricard/10x_gastrulation_DNMTs/utils.R")
}

# I/O

# Options
opts$remove.ExE.celltypes <- FALSE
opts$min.cells <- 10

## START TEST ##
args$sce <- io$sce
args$metadata <- io$metadata
# args$samples <- opts$samples[1]
# args$features <- 1000
# args$npcs <- 30
# args$test <- TRUE
# args$colour_by <- c("celltype.mapped")
# args$vars.to.regress <- c("nFeature_RNA","percent.mt")
# args$batch.correction <- c("sample")
# args$outdir <- paste0(io$basedir,"/results/dimensionality_reduction/test")
## END TEST ##

# if (isTRUE(args$test)) print("Test mode activated...")

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_QC==TRUE & sample%in%args$samples]

if (opts$remove.ExE.celltypes) {
  sample_metadata <- sample_metadata %>%
    .[!celltype.mapped%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

sample_metadata <- sample_metadata %>%
  .[,N:=.N,by=c("celltype.mapped")] %>% .[N>opts$min.cells] %>% .[,N:=NULL]

table(sample_metadata$sample)
table(sample_metadata$celltype.mapped)

###################
## Sanity checks ##
###################

stopifnot(args$colour_by %in% colnames(sample_metadata))

stopifnot(unique(sample_metadata$celltype.mapped) %in% names(opts$celltype.colors))

if (length(args$batch.correction)>0) {
  stopifnot(args$batch.correction%in%colnames(sample_metadata))
  if (length(unique(sample_metadata[[args$batch.correction]]))==1) {
    message(sprintf("There is a single level for %s, no batch correction applied",args$batch.correction))
    args$batch.correction <- NULL
  } else {
    library(batchelor)
  }
}

if (length(args$vars.to.regress)>0) {
  stopifnot(args$vars.to.regress%in%colnames(sample_metadata))
}


###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- load_SingleCellExperiment(args$sce, cells=sample_metadata$cell, normalise = TRUE)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

#######################
## Feature selection ##
#######################

if (length(args$batch.correction)>0) {
  decomp <- modelGeneVar(sce, block=colData(sce)[[args$batch.correction]])
} else {
  decomp <- modelGeneVar(sce)
  
}
decomp <- decomp[decomp$mean > 0.01,]
hvgs <- decomp[order(decomp$FDR),] %>% head(n=args$features) %>% rownames

# Subset SingleCellExperiment
sce_filt <- sce[hvgs,]

############################
## Regress out covariates ##
############################

if (length(args$vars.to.regress)>0) {
  print(sprintf("Regressing out variables: %s", paste(args$vars.to.regress,collapse=" ")))
  logcounts(sce_filt) <- RegressOutMatrix(
    mtx = logcounts(sce_filt), 
    covariates = colData(sce_filt)[,args$vars.to.regress,drop=F]
  )
}

############################
## PCA + Batch correction ##
############################

if (length(args$batch.correction)>0) {
  print(sprintf("Applying MNN batch correction for variable: %s", args$batch.correction))
  outfile <- sprintf("%s/%s_pca_features%d_pcs%d_batchcorrectionby%s.txt.gz",args$outdir, paste(args$samples,collapse="-"), args$features, args$npcs,paste(args$batch.correction,collapse="-"))
  pca <- multiBatchPCA(sce_filt, batch = colData(sce_filt)[[args$batch.correction]], d = args$npcs)
  pca.corrected <- reducedMNN(pca)$corrected
  colnames(pca.corrected) <- paste0("PC",1:ncol(pca.corrected))
  reducedDim(sce_filt, "PCA") <- pca.corrected
} else {
  outfile <- sprintf("%s/%s_pca_features%d_pcs%d.txt.gz",args$outdir, paste(args$samples,collapse="-"), args$features, args$npcs)
  sce_filt <- runPCA(sce_filt, ncomponents = args$npcs, ntop=args$features)  
}

# Save PCA coordinates
pca.dt <- reducedDim(sce_filt,"PCA") %>% round(3) %>% as.data.table(keep.rownames = T) %>% setnames("rn","cell")
fwrite(pca.dt, outfile)

##########
## UMAP ##
##########

for (i in args$n_neighbors) {
  for (j in args$min_dist) {
    
    # Run
    set.seed(args$seed)
    sce_filt <- runUMAP(sce_filt, dimred="PCA", n_neighbors = i, min_dist = j)
    
    # Fetch UMAP coordinates
    umap.dt <- reducedDim(sce_filt,"UMAP") %>% as.data.table %>% 
      .[,cell:=colnames(sce_filt)] %>%
      setnames(c("UMAP1","UMAP2","cell"))
    
    # Plot
    to.plot <- reducedDim(sce_filt,"UMAP") %>% as.data.table %>% 
      .[,cell:=colnames(sce_filt)] %>%
      merge(sample_metadata, by="cell")
    
    for (k in args$colour_by) {
      
      p <- ggplot(to.plot, aes_string(x="V1", y="V2", fill=k)) +
        geom_point(size=1.5, shape=21, stroke=0.05) +
        theme_classic() +
        theme(
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()
        )
      
      if (k=="celltype.mapped") {
        p <- p + scale_fill_manual(values=opts$celltype.colors) +
          theme(
            legend.position="none",
            legend.title=element_blank()
          )
      }
      
      # Save UMAP plot
      outfile <- sprintf("%s/%s_umap_features%d_pcs%d_neigh%d_dist%s_%s.pdf",args$outdir, paste(args$samples,collapse="-"), args$features, args$npcs, i, j, k)
      pdf(outfile, width=7, height=5)
      print(p)
      dev.off()

    }
    
    # Save UMAP coordinates
    outfile <- sprintf("%s/%s_umap_features%d_pcs%d_neigh%d_dist%s.txt.gz",args$outdir, paste(args$samples,collapse="-"), args$features, args$npcs, i, j)
    fwrite(umap.dt, outfile)
  }
}

