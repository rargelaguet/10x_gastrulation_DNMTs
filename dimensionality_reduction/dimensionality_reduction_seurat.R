here::i_am("dimensionality_reduction/dimensionality_reduction_seurat.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(Seurat))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--seurat',             type="character",                               help='Seurat object file')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--classes',       type="character",  default="all",  nargs='+',  help='Classes to plot')
p$add_argument('--features',        type="integer",    default=1000,                help='Number of features')
p$add_argument('--npcs',            type="integer",    default=30,                  help='Number of PCs')
p$add_argument('--vars_to_regress', type="character",                nargs='+',     help='Metadata columns to regress out')
p$add_argument('--n_neighbors',     type="integer",    default=30,     help='(UMAP) Number of neighbours')
p$add_argument('--min_dist',        type="double",     default=0.3,     help='(UMAP) Minimum distance')
p$add_argument('--colour_by',       type="character",  default="celltype.mapped",  nargs='+',  help='Metadata columns to colour the UMAP by')
p$add_argument('--seed',            type="integer",    default=42,                  help='Random seed')
p$add_argument('--outdir',          type="character",                               help='Output file')
p$add_argument('--remove_ExE_cells', action="store_true",                                 help='Remove ExE cells?')
p$add_argument('--batch_correction', type="character",                               help='Metadata column to apply batch correction on')
# p$add_argument('--SCTransform', action="store_true",                                 help='Remove ExE cells?')

args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
# args$seurat <- file.path(io$basedir,"processed_all/seurat.rds")
# args$metadata <- file.path(io$basedir,"results_all/mapping/sample_metadata_after_mapping.txt.gz")
# args$classes <- "Dnmt1_KO"
# args$features <- 2500
# args$npcs <- 30
# args$colour_by <- c("celltype.mapped","sample","nFeature_RNA","dataset")
# args$batch_correction <- "dataset"
# args$vars_to_regress <- c("nFeature_RNA","mit_percent_RNA")
# args$n_neighbors <- 25
# args$min_dist <- 0.5
# args$seed <- 42
# args$outdir <- file.path(io$basedir,"results_all/dimensionality_reduction/seurat")
# args$remove_ExE_cells <- TRUE
## END TEST ##

# if (isTRUE(args$test)) print("Test mode activated...")

# I/O
dir.create(args$outdir, showWarnings = F)

# Options
if (args$classes[1]=="all") {
  args$classes <- opts$classes
} else {
  stopifnot(args$classes%in%opts$classes)
}

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  .[,dataset:=ifelse(grepl("Grosswendt",sample),"Grosswendt","This data set")] %>%
  .[pass_rnaQC==TRUE & class%in%args$classes]

if (args$remove_ExE_cells) {
  print("Removing ExE cells...")
  sample_metadata <- sample_metadata %>%
    .[!celltype.mapped%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

table(sample_metadata$class)
table(sample_metadata$dataset)
table(sample_metadata$celltype.mapped)

###################
## Sanity checks ##
###################

stopifnot(args$colour_by %in% colnames(sample_metadata))
# stopifnot(unique(sample_metadata$celltype.mapped) %in% names(opts$celltype.colors))

if (length(args$batch_correction)>0) {
  stopifnot(args$batch_correction%in%colnames(sample_metadata))
  if (length(unique(sample_metadata[[args$batch_correction]]))==1) {
    message(sprintf("There is a single level for %s, no batch correction applied",args$batch_correction))
    args$batch_correction <- NULL
  }
}

if (length(args$vars_to_regress)>0) {
  stopifnot(args$vars_to_regress%in%colnames(sample_metadata))
}


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

#################
## SCTransform ##
#################

#   variable.features.n = args$features, 
#   # vars_to_regress = c("nFeature_RNA","mitochondrial_percent_RNA"),
#   vars_to_regress = args$vars_to_regress,
#   verbose = FALSE
# )

#######################
## Feature selection ##
#######################

if (is.null(args$batch_correction)) {
  seurat <- FindVariableFeatures(seurat, selection.method = 'vst', nfeatures = args$features)
  # seurat <- FindVariableFeatures(seurat, selection.method = 'vst', nfeatures = 1000, assay = "SCT")
}

###########################################
## Scale data and regress out covariates ##
###########################################

if (is.null(args$batch_correction)) {
# seurat <- ScaleData(seurat, features=var.genes, vars_to_regress=c("nFeature_RNA","mitochondrial_percent_RNA"))
  seurat <- ScaleData(seurat, 
    features = VariableFeatures(seurat), 
    vars.to.regress = args$vars_to_regress, 
    verbose = FALSE
  )
}


######################
## Batch correction ##
######################

if (length(args$batch_correction)>0) {
  seurat.list <- SplitObject(seurat, split.by = args$batch_correction)

  for (i in 1:length(seurat.list)) {
    
    # Normalisation
    seurat.list[[i]] <- NormalizeData(
      object = seurat.list[[i]], 
      verbose = FALSE
    )
    
    # Feature selection
    seurat.list[[i]] <- FindVariableFeatures(
      object = seurat.list[[i]],
      selection.method = "vst",
      nfeatures = 2000, 
      verbose = FALSE
    )
  }

  anchors <- FindIntegrationAnchors(object.list = seurat.list, dims = 1:30)

  seurat <- IntegrateData(anchorset = anchors, dims = 1:30)

  DefaultAssay(seurat) <- "integrated"
  seurat <- ScaleData(seurat, verbose = FALSE)
  
  rm(anchors,seurat.list)
}

#########
## PCA ##
#########

seurat <- RunPCA(seurat, features = VariableFeatures(seurat), npcs = args$npcs, verbose = FALSE)

# Save PCA coordinates
pca.dt <- seurat@reductions[["pca"]]@cell.embeddings %>% round(3) %>% as.data.table(keep.rownames = T) %>% setnames("rn","cell")
fwrite(pca.dt, file.path(args$outdir,sprintf("pca_features%d_pcs%d.txt.gz", args$features, args$npcs)))

##########
## UMAP ##
##########

# Run
set.seed(args$seed)
seurat <- RunUMAP(seurat, 
  dims = 1:args$npcs,
  reduction = "pca",
  n.neighbors = args$n_neighbors,
  min.dist = args$min_dist
)

# Fetch UMAP coordinates
umap.dt <- seurat@reductions[["umap"]]@cell.embeddings %>% round(3) %>% as.data.table %>% 
  .[,cell:=colnames(seurat)] %>%
  setnames(c("UMAP1","UMAP2","cell")) %>% 
  .[,c("cell","UMAP1","UMAP2")]

# Save UMAP coordinates
fwrite(umap.dt, file.path(args$outdir,sprintf("umap_features%d_pcs%d_neigh%d_dist%s.txt.gz",args$features, args$npcs, args$n_neighbors, args$min_dist)))

##########
## Plot ##
##########

pt.size <- ifelse(ncol(seurat)>=1e4,1,1.25)

for (i in args$colour_by) {
  
  to.plot <- umap.dt %>%
    setnames(c("cell","V1","V2")) %>%
    merge(sample_metadata, by="cell")
  
  if (is.numeric(to.plot[[i]])) {
    if ((max(to.plot[[i]],na.rm=T) - min(to.plot[[i]],na.rm=T)) > 1000) {
      to.plot[[i]] <- log10(to.plot[[i]]+1)
      to.plot %>% setnames(i,paste0(i,"_log10")); i <- paste0(i,"_log10")
    }
  }
  
  p <- ggplot(to.plot, aes_string(x="V1", y="V2", fill=i)) +
    geom_point(size=pt.size, shape=21, stroke=0.05) +
    theme_classic() +
    ggplot_theme_NoAxes()
  

  if (grepl("celltype",i)) {
    p <- p + scale_fill_manual(values=opts$celltype.colors) + NoLegend()
  }
  
  if (grepl("stage",i)) {
    p <- p + scale_fill_manual(values=opts$stage.colors) + NoLegend()
  }  
  
  if (grepl("sample",i)) {
    p <- p + theme(
      legend.position = "none",
      legend.title = element_blank()
    )
  }  
  
  # Define colormap
  if (is.numeric(to.plot[[i]])) {
    p <- p + scale_fill_gradientn(colours = terrain.colors(10))
  }
  

  # Save UMAP plot
  outfile <- file.path(args$outdir,sprintf("umap_features%d_pcs%d_neigh%d_dist%s_%s.pdf", args$features, args$npcs, args$n_neighbors, args$min_dist, i))
  pdf(outfile, width=7, height=5)
  print(p)
  dev.off()
}


# Completion token
file.create(file.path(args$outdir,"completed.txt"))