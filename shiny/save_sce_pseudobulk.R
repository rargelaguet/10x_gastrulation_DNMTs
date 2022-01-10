here::i_am("plot_individual_genes/pseudobulk/plot_individual_genes_pseudobulk.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$outfile <- "/Users/argelagr/data/shiny_dnmt_tet/SingleCellExperiment_pseudobulk_class_celltype_dataset.rds"
io$sce.pseudobulk <- file.path(io$basedir,"results_all/pseudobulk/SingleCellExperiment_pseudobulk_class_celltype_dataset.rds")

# Define options
opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  # "PGC",
  "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  "Mixed_mesoderm",
  "Intermediate_mesoderm",
  "Caudal_Mesoderm",
  "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "Pharyngeal_mesoderm",
  "Cardiomyocytes",
  "Allantois",
  "ExE_mesoderm",
  "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  # "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm"
  # "ExE_ectoderm",
  # "Parietal_endoderm"
)

opts$classes <- c(
  # "Dnmt3a_HET_Dnmt3b_KO",
  # "Dnmt3a_HET_Dnmt3b_WT",
  # "Dnmt3a_KO_Dnmt3b_HET",
  "WT",
  "Dnmt3a_KO",
  "Dnmt3b_KO",
  "Dnmt1_KO",
  "Dnmt3ab_KO"
)

# opts$rename_celltypes <- c(
#   "Erythroid3" = "Erythroid",
#   "Erythroid2" = "Erythroid",
#   "Erythroid1" = "Erythroid",
#   "Blood_progenitors_1" = "Blood_progenitors",
#   "Blood_progenitors_2" = "Blood_progenitors",
#   "Anterior_Primitive_Streak" = "Primitive_Streak"
# )

opts$min.cells <- 30

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- readRDS(io$sce.pseudobulk)

# Add metadata
sce$class <- stringr::str_split(colnames(sce), pattern = "-") %>% map_chr(1)
sce$celltype <- stringr::str_split(colnames(sce), pattern = "-") %>% map_chr(2)
sce$dataset <- stringr::str_split(colnames(sce), pattern = "-") %>% map_chr(3)

# Filter by minimum number of cells
sce <- sce[,names(which(metadata(sce)$n_cells>=opts$min.cells))]

# Filter
sce <- sce[,sce$celltype%in%opts$celltypes]
sce <- sce[,sce$class%in%opts$classes]

genes <- fread("/Users/argelagr/data/shiny_dnmt_tet/genes.txt", header = F)[[1]]
sce <- sce[genes,]

##########
## Save ##
##########

assays(sce)["counts"] <- NULL
saveRDS(sce, io$outfile)
