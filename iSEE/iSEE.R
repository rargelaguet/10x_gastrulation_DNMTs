library(SingleCellExperiment)
library(scran)
library(scater)
library(shiny)
library(iSEE)

#####################
## Define settings ##
#####################

# Load default setings

source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")

# Define I/O
# io$outdir <- paste0(io$basedir,"/results/dimensionality_reduction")

# Define options
opts$batches <- c(
  # E12.5  
  # "E125_DNMT3A_HET_A_L001",
  # "E125_DNMT3A_HET_A_L003",
  # "E125_DNMT3A_KO_B_L002",
  # "E125_DNMT3A_KO_E_L004",
  # "A_E12_5_D3a_Het_L001",
  # "B_E12_5_D3a_KO_L002",
  
  
  # E8.5  
  # "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001",
  # "SIGAB6_E85_3_Dnmt3aWT_Dnmt3b_WT_L002",
  # "SIGAC6_E85_5_Dnmt3aKO_Dnmt3b_Het_L003",
  # "SIGAD6_E85_8_Dnmt3aHet_Dnmt3b_KO_L004",
  "15_E8_5_D3A_WT_D3B_WT_L007"
  # "17_E8_5_D3A_KO_D3B_WT_L008",
  # "2_E8_5_D3A_WT_D3B_KO_L003",
  # "3_E8_5_D3A_HET_D3B_WT_L004",
  # "7_E8_5_D3A_WT_D3B_KO_L005",
  # "8_E8_5_D3A_KO_D3B_KO_L006",
  # "E8_5_Dnmt1_KO_male_SIGAC8_L001",
  # "E8_5_Dnmt1_KO_male_SIGAD8_L002",
  # "E8_5_Dnmt1_KO_male_SIGAE8_L003",
  # "E8_5_Dnmt1_WT_female_SIGAB8_L004",
  # "E8_5_Dnmt1_WT_female_SIGAF8_L005",
  # "E8_5_Dnmt3ab_WT_female_SIGAA8_L006",
  # "SIGAH10_Dnmt3ab_WT_L002",
  # "SIGAH11_Dnmt3ab_WT_L003",
  # "SIGAH9_Dnmt3a_KO_Dnmt3b_Het_L001"
)


############################
## Update sample metadata ##
############################

sample_metadata <- fread(io$metadata) %>%
  .[pass_QC==TRUE & batch%in%opts$batches]
table(sample_metadata$batch)

###############
## Load data ##
###############

# Load RNA expression data as SingleCellExperiment object
sce <- readRDS(io$sce)[,sample_metadata$cell]

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

# Add gene metadata as rowData
# rowData(sce)

#######################
## Feature selection ##
#######################

decomp <- modelGeneVar(sce)
decomp <- decomp[decomp$mean > 0.01,]
hvgs <- rownames(decomp)[decomp$p.value <= 0.05]

# Subset SingleCellExperiment
sce_filt <- sce[hvgs,]
dim(sce_filt)

##############################
## Dimensionality reduction ##
##############################

data <- scale(t(logcounts(sce_filt)), center = T, scale = F)

# PCA
npcs <- 30
reducedDim(sce_filt, "PCA") <- irlba::prcomp_irlba(data, n=npcs)$x#[,1:npcs]

# Filter PCA solution
# reducedDim(sce_filt, "PCA") <- reducedDim(sce_filt, "PCA")[,-3]

# UMAP
set.seed(42)
sce_filt <- runUMAP(sce_filt, dimred="PCA", n_neighbors = 30, min_dist = 0.3)

#######################
## Define color maps ##
#######################

celltype.colors <- opts$celltype.colors[names(opts$celltype.colors)%in%unique(sce_filt$celltype.mapped)]
all(unique(sce_filt$celltype.mapped) %in% names(celltype.colors))
all(names(celltype.colors)%in%unique(sce_filt$celltype.mapped))
sce_filt$celltype.mapped <- factor(sce_filt$celltype.mapped, levels=names(celltype.colors))
celltype_color_fun <- function(){
  return(celltype.colors)
}

categorical_color_fun <- function(n){
  return(RColorBrewer::brewer.pal(n, "Set2"))
}

# Define color maps
ecm <- ExperimentColorMap(
  # List of colormaps for assays.
  # assays = list(
  #   counts = viridis::viridis,
  #   cufflinks_fpkm = fpkm_color_fun
  # ),
  colData = list(
    celltype.mapped = celltype_color_fun
  ),
  # Colormaps applied to all undefined continuous assays
  all_continuous = list(
    assays = viridis::viridis
  ),
  # Colormaps applied to all undefined categorical assays
  all_discrete = list(
    assays = categorical_color_fun
  )
  # Colormap applied to all undefined categorical covariates.
  # global_discrete <- list()
  # Colormap applied to all undefined continuous covariates.
  # global_continuous <- list()
)

##############
## Run iSEE ##
##############

app <- iSEE(sce_filt, colormap = ecm)
runApp(app)
