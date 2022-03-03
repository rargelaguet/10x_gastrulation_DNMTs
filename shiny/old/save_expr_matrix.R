here::i_am("shiny/save_expr_matrix.R")

library(HDF5Array)

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O ##
io$sce <- file.path(io$basedir,"processed/SingleCellExperiment.rds")
io$outdir <- "/Users/argelagr/data/shiny_dnmt_tet"

## Define options ##

# Define cell types to plot
opts$celltypes = c(
	"Epiblast",
	"Primitive_Streak",
	"Caudal_epiblast",
	"PGC",
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
	# "Blood_progenitors",
	"Blood_progenitors_1",
	"Blood_progenitors_2",
	# "Erythroid",
	"Erythroid1",
	"Erythroid2",
	"Erythroid3",
	"NMP",
	"Rostral_neurectoderm",
	"Caudal_neurectoderm",
	"Neural_crest",
	"Forebrain_Midbrain_Hindbrain",
	"Spinal_cord",
	"Surface_ectoderm",
	"Visceral_endoderm",
	"ExE_endoderm",
	"ExE_ectoderm",
	"Parietal_endoderm"
)

# Define classes to plot
opts$classes <- c(
  "WT",
  "Dnmt3a_KO",
  "Dnmt3b_KO",
  "Dnmt1_KO"
  # "Dnmt3ab_KO"
)

# opts$rename_celltypes <- c(
#   "Erythroid3" = "Erythroid",
#   "Erythroid2" = "Erythroid",
#   "Erythroid1" = "Erythroid",
#   "Blood_progenitors_1" = "Blood_progenitors",
#   "Blood_progenitors_2" = "Blood_progenitors",
#   "Anterior_Primitive_Streak" = "Primitive_Streak"
# )

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  .[,dataset:=ifelse(grepl("Grosswendt",sample),"CRISPR","KO")] %>%
  # .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,opts$rename_celltypes)] %>%
  .[pass_rnaQC==TRUE & celltype.mapped%in%opts$celltypes & class%in%opts$classes]

table(sample_metadata$dataset)
table(sample_metadata$class)
table(sample_metadata$celltype.mapped)

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- load_SingleCellExperiment(io$sce, cells=sample_metadata$cell, normalise = TRUE)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

#################
## Filter data ##
#################

sce <- sce[!grepl("*Rik|^Gm|^mt-|^Rps|^Rpl|^Olf",rownames(sce)),]

# Remove genes that are not expressed
sce <- sce[sparseMatrixStats::rowVars(logcounts(sce))>0,]

##########
## Save ##
##########

# Save metadata
cols <- c(c("cell","alias","class","dataset","stage","nFeature_RNA", "mit_percent_RNA", "rib_percent_RNA", "celltype.mapped", "celltype.score", "closest.cell"))
sample_metadata_to_save <- sample_metadata[,..cols]
fwrite(sample_metadata_to_save, file.path(io$outdir,"cell_metadata.txt.gz"), quote=F, sep="\t", na="NA")

# Save gene names
write.table(rownames(sce), paste0(io$outdir,"/genes.txt"), row.names = F, col.names = F, quote=F)

# Save cell names
write.table(colnames(sce), paste0(io$outdir,"/cells_rna.txt"), row.names = F, col.names = F, quote=F)

# Save expression matrix
io$hdf5.outfile <- file.path(io$outdir,"rna_expr.hdf5")
if(file.exists(io$hdf5.outfile)) { file.remove(io$hdf5.outfile) }
writeHDF5Array(x = DelayedArray(round(logcounts(sce),2)), file = io$hdf5.outfile, name = "expr_logcounts", verbose = TRUE)
