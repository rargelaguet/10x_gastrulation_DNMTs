here::i_am("trajectories/infer_trajectory.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(destiny))
suppressPackageStartupMessages(library(reticulate))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--metadata',  type="character",              help='') 
p$add_argument('--sce',  type="character",              help='') 
p$add_argument('--trajectory_name',  type="character",              help='') 
p$add_argument('--celltype_label',  type="character",              help='') 
# p$add_argument('--genes_to_plot',  type="character", nargs="+", help='') 
p$add_argument('--outdir',          type="character",                help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
args <- list()
args$sce <-file.path(io$basedir,"processed_new/SingleCellExperiment.rds")
args$metadata <- file.path(io$basedir,"results_new/mapping/sample_metadata_after_mapping.txt.gz")
args$trajectory_name <- "NMP"
args$celltype_label <- "celltype.mapped"
args$vars_to_regress <- c("nFeature_RNA")
args$outdir <- file.path(io$basedir,"results_new/trajectories/NMP")
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings=F)

# Options
opts$celltype_trajectory_dic <- list(
  "blood" = c("Haematoendothelial_progenitors", "Blood_progenitors_1", "Blood_progenitors_2", "Erythroid1", "Erythroid2", "Erythroid3"),
  "ectoderm" = c("Epiblast", "Rostral_neurectoderm", "Forebrain_Midbrain_Hindbrain"),
  "endoderm" = c("Epiblast", "Anterior_Primitive_Streak", "Def._endoderm", "Gut"),
  "mesoderm" = c("Epiblast", "Primitive_Streak", "Nascent_mesoderm"),
  # "NMP" = c("Caudal_epiblast", "Somitic_mesoderm", "Spinal_cord","Caudal_Mesoderm","NMP")
  "NMP" = c("Epiblast","Primitive_Streak","Caudal_epiblast","Caudal_Mesoderm","NMP")
)

stopifnot(args$trajectory_name%in%names(opts$celltype_trajectory_dic))
opts$celltypes <- opts$celltype_trajectory_dic[[args$trajectory_name]]

# opts$genes2plot <- list(
#   "blood" = c("Hbb-y","Etv2"),
#   "ectoderm" = c("Utf1","Crabp1"),
#   "endoderm" = c("Pou5f1","Krt8"),
#   "mesoderm" = c("Dnmt3b","Mesp1")
# )
# stopifnot(args$trajectory_name%in%names(opts$genes2plot))
# opts$genes_to_plot <- opts$genes2plot[[args$trajectory_name]]
opts$genes_to_plot <- grep("Hox",unique(fread(io$atlas.marker_genes)[,gene]),value=T)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_rnaQC==TRUE & class%in%c("E8.5_WT") & nFeature_RNA>=2500]

stopifnot(args$celltype_label%in%colnames(sample_metadata))
sample_metadata <- sample_metadata %>%
  .[,celltype:=eval(as.name(args$celltype_label))] %>%
  .[celltype%in%opts$celltypes] %>%
  .[,celltype:=factor(celltype,levels=opts$celltypes)]

table(sample_metadata$celltype)

# Save
fwrite(sample_metadata, file.path(args$outdir,sprintf("%s_sample_metadata.txt.gz",args$trajectory_name)), sep="\t", na="NA")

#########################
## Load RNA expression ##
#########################
  
sce <- load_SingleCellExperiment(
  file = args$sce, 
  normalise = TRUE, 
  cells = sample_metadata$cell, 
  remove_non_expressed_genes = TRUE
)

colData(sce) <- sample_metadata %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
  .[colnames(sce),] %>% DataFrame()

#######################
## Feature selection ##
#######################

decomp <- modelGeneVar(sce)
decomp <- decomp[decomp$mean > 0.01,]
hvgs <- rownames(decomp)[decomp$p.value <= 0.01]

# Subset HVGs
sce_filt <- sce[hvgs,]
dim(sce_filt)

#####################################
## Regress out technical variables ##
#####################################

if (length(args$vars_to_regress)>0) {
  print(sprintf("Regressing out variables: %s", paste(args$vars_to_regress,collapse=" ")))
  logcounts(sce_filt) <- RegressOutMatrix(
    mtx = logcounts(sce_filt),
    covariates = colData(sce_filt)[,args$vars_to_regress,drop=F]
  )
}


#########
## PCA ##
#########

sce_filt <- runPCA(sce_filt, ncomponents = 5, ntop=nrow(sce_filt))

# Plot PCA
pdf(file.path(args$outdir,sprintf("%s_pca_celltype.pdf",args$trajectory_name)), width=8, height=5)
plotPCA(sce_filt, colour_by="celltype", ncomponents = c(1,2)) +
  scale_colour_manual(values=opts$celltype.colors) + theme(legend.position="none")
dev.off()

######################
## Batch correction ##
######################

# args$batch_correction <- "stage"
# if (length(args$batch_correction)>0) {
#   suppressPackageStartupMessages(library(batchelor))
#   print(sprintf("Applying MNN batch correction for variable: %s", args$batch_correction))
#   pca <- multiBatchPCA(sce_filt, batch = colData(sce_filt)[[args$batch_correction]], d = args$npcs)
#   pca.corrected <- reducedMNN(pca)$corrected
#   colnames(pca.corrected) <- paste0("PC",1:ncol(pca.corrected))
#   reducedDim(sce_filt, "PCA") <- pca.corrected
# }

###################
## Diffusion map ##
###################

set.seed(42)
dm <- DiffusionMap(sce_filt, n_pcs=2)

# Add to the SingleCellExperiment object
reducedDim(sce_filt, "DiffusionMap") <- dm@eigenvectors[,c(1,2)]

# Plot
pdf(file.path(args$outdir,sprintf("%s_diffmap_celltype.pdf",args$trajectory_name)), width=8, height=5)
plotReducedDim(sce_filt, dimred = "DiffusionMap", colour_by="celltype", ncomponents = c(1,2)) +
  scale_colour_manual(values=opts$celltype.colors) + theme(legend.position="none")
dev.off()


##########
## UMAP ##
##########

set.seed(42)

sce_filt <- runUMAP(sce_filt, dimred="PCA", n_neighbors = 15, min_dist = 0.3)

# Plot
pdf(file.path(args$outdir,sprintf("%s_umap_celltype.pdf",args$trajectory_name)), width=8, height=5)
plotReducedDim(sce_filt, dimred = "UMAP", colour_by="celltype", ncomponents = c(1,2)) +
  scale_colour_manual(values=opts$celltype.colors) + theme(legend.position="none")
dev.off()

####################
## Prepare output ##
####################

pseudotime.dt <- data.table(
  cell = colnames(sce_filt),
  DC1 = reducedDim(sce_filt,"DiffusionMap")[,"DC1"],
  DC2 = reducedDim(sce_filt,"DiffusionMap")[,"DC2"]
)

# Save
fwrite(pseudotime.dt, file.path(args$outdir,sprintf("/%s_trajectory.txt.gz",args$trajectory_name)))

#########################################################################
## Scatterplots of pseudotime values versus expression of marker genes ##
#########################################################################

# Denoise
# pca.rna <- fread(io$pca.rna) %>% matrix.please %>% .[colnames(sce),]
# assay(sce_filt,"logcounts_denoised") <- smoother_aggregate_nearest_nb(mat=as.matrix(logcounts(sce_filt)), D=pdist(pca.rna), k=25)

rna.dt <- data.table(as.matrix(logcounts(sce)[opts$genes_to_plot,]), keep.rownames = T) %>%
  setnames("rn","gene") %>%
  melt(id.vars="gene", variable.name="cell", value.name="expr") %>%
  merge(sample_metadata[,c("cell","celltype.mapped")])

opts$genes_to_plot <- c("Hoxc9","Hoxb9","Hoxa9","Hoxd9", "Hoxc8","Hoxb8", "Hoxa7","Hoxb7", "Hoxb6","Hoxc6")

to.plot <- pseudotime.dt %>%
  merge(rna.dt[gene%in%opts$genes_to_plot],by="cell") %>%
  .[,gene:=factor(gene,levels=opts$genes_to_plot)] %>%
  .[,expr:=minmax.normalisation(expr),by="gene"]

p <- ggplot(to.plot, aes_string(x="DC1", y="DC2")) +
  geom_point(aes(fill=expr), size=1, shape=21, stroke=0.1, alpha=0.75) +
  # viridis::scale_fill_viridis() +
  scale_fill_gradient2(low = "gray50", mid="gray90", high = "red") +
  facet_wrap(~gene, scales="free_y") +
  # labs(x=sprintf("Pseudotime (%s)",i), y="Gene expression") +
  theme_classic() +
  guides(fill="none", color="none") +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    # axis.text.y = element_text(color="black"),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_text(size=rel(1), color="black"),
    legend.position="right",
    legend.title = element_blank()
  )

pdf(file.path(args$outdir,sprintf("%s_gene_expression.pdf",args$trajectory_name)), width=12, height=8)
print(p)
dev.off()


###############################
## Save SingleCellExperiment ##
###############################

logcounts(sce) <- NULL
saveRDS(sce, file.path(args$outdir,sprintf("%s_SingleCellExperiment.rds",args$trajectory_name)))

