here::i_am("shiny/save_sce_pseudobulk.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$outfile <- "/Users/argelagr/data/shiny_dnmt_tet/SingleCellExperiment_pseudobulk_class_sample_celltype.rds"
io$sce.pseudobulk <- file.path(io$basedir,"results/pseudobulk/SingleCellExperiment_pseudobulk_class_sample_celltype.rds")

# Define options
opts$celltypes <- c(
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
  "WT",
  "Dnmt3a_KO",
  "Dnmt3b_KO",
  "Dnmt1_KO"
  # "Dnmt3ab_KO"
)


opts$min.cells <- 30

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- readRDS(io$sce.pseudobulk)

# Add metadata
sce$class <- stringr::str_split(colnames(sce), pattern = "-") %>% map_chr(1)
sce$sample <- stringr::str_split(colnames(sce), pattern = "-") %>% map_chr(2)
sce$celltype <- stringr::str_split(colnames(sce), pattern = "-") %>% map_chr(3)

# Add data set info
sce$dataset <- ifelse(grepl("Grosswendt",colnames(sce)),"CRISPR","KO")

# Filter by minimum number of cells
metadata(sce)$n_cells <- metadata(sce)$n_cells[metadata(sce)$n_cells>=opts$min.cells]
sce <- sce[,names(metadata(sce)$n_cells)]

# Filter classes and celltypes manually
sce <- sce[,sce$celltype%in%opts$celltypes]
sce <- sce[,sce$class%in%opts$classes]
metadata(sce)$n_cells <- metadata(sce)$n_cells[stringr::str_split(names(metadata(sce)$n_cells), pattern = "-") %>% map_chr(3) %in% opts$celltypes]
metadata(sce)$n_cells <- metadata(sce)$n_cells[stringr::str_split(names(metadata(sce)$n_cells), pattern = "-") %>% map_chr(1) %in% opts$classes]

# Filter celltypes with not enough observations
stats.dt <- data.table(class=sce$class, celltype=sce$celltype, sample=sce$sample, id=colnames(sce)) %>% 
  .[,N:=length(unique(class)),by=c("celltype")] %>% .[N>=3]
sce <- sce[,stats.dt$id]

##################
## Filter genes ##
##################

genes <- fread("/Users/argelagr/data/shiny_dnmt_tet/genes.txt", header = F)[[1]]
sce <- sce[genes,]

##########
## Save ##
##########

assays(sce)["counts"] <- NULL
saveRDS(sce, io$outfile)
