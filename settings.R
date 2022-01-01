suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(argparse))

#########
## I/O ##
#########

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/10x_gastrulation_DNMTs"
  io$atlas.basedir <- "/Users/ricard/data/gastrulation10x"
  io$gene_metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/all_genes/Mmusculus_genes_BioMart.87.txt"
} else if (grepl("ebi",Sys.info()['nodename'])) {
  io$basedir <- "/hps/nobackup2/research/stegle/users/ricard/10x_gastrulation_DNMTs"
  io$atlas.basedir <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x"
  io$gene_metadata <- "/hps/nobackup2/research/stegle/users/ricard/ensembl/mouse/v87/BioMart/all_genes/Mmusculus_genes_BioMart.87.txt"
} else if (Sys.info()[['nodename']]=="BI2404M") {
  io$basedir <- "/Users/argelagr/data/10x_gastrulation_DNMTs"
  io$atlas.basedir <- "/Users/argelagr/data/gastrulation10x"
  io$gene_metadata <- "/Users/argelagr/data/ensembl/mouse/v87/BioMart/all_genes/Mmusculus_genes_BioMart.87.txt"
} else if (grepl("pebble|headstone", Sys.info()['nodename'])) {
  if (grepl("argelag", Sys.info()['effective_user'])) {
    io$basedir <- "/bi/group/reik/ricard/data/10x_gastrulation_DNMTs"
    io$atlas.basedir <- "/bi/group/reik/ricard/data/pijuansala2019_gastrulation10x"
    io$gene_metadata <- "/bi/group/reik/ricard/data/ensembl/mouse/v87/BioMart/all_genes/Mmusculus_genes_BioMart.87.txt"
  }
} else {
  stop("Computer not recognised")
}

io$metadata <- paste0(io$basedir,"/sample_metadata.txt.gz")
io$seurat <- paste0(io$basedir,"/processed_new/seurat.rds")
io$sce <- paste0(io$basedir,"/processed_new/SingleCellExperiment.rds")
# io$sex <- paste0(io$basedir,"/results/sex/sex_assignment.txt.gz")

# Atlas information
io$atlas.metadata <- paste0(io$atlas.basedir,"/sample_metadata.txt.gz")
io$atlas.marker_genes <- paste0(io$atlas.basedir,"/results/marker_genes/all_stages/marker_genes.txt.gz")
io$atlas.differential <- paste0(io$atlas.basedir,"/results/differential")
# io$atlas.average_expression_per_celltype <- paste0(io$atlas.basedir,"/results/marker_genes/avg_expr_per_celltype_and_gene.txt.gz")
io$atlas.sce <- paste0(io$atlas.basedir,"/processed/SingleCellExperiment.rds")


#############
## Options ##
#############

opts <- list()

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
	"Blood_progenitors_1",
	"Blood_progenitors_2",
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

opts$celltype.colors = c(
  "Epiblast" = "#635547",
  "Primitive_Streak" = "#DABE99",
  "Caudal_epiblast" = "#9e6762",
  "PGC" = "#FACB12",
  "Anterior_Primitive_Streak" = "#c19f70",
  "Notochord" = "#0F4A9C",
  "Def._endoderm" = "#F397C0",
  "Gut" = "#EF5A9D",
  "Nascent_mesoderm" = "#C594BF",
  "Mixed_mesoderm" = "#DFCDE4",
  "Intermediate_mesoderm" = "#139992",
  "Caudal_Mesoderm" = "#3F84AA",
  "Paraxial_mesoderm" = "#8DB5CE",
  "Somitic_mesoderm" = "#005579",
  "Pharyngeal_mesoderm" = "#C9EBFB",
  "Cardiomyocytes" = "#B51D8D",
  "Allantois" = "#532C8A",
  "ExE_mesoderm" = "#8870ad",
  "Mesenchyme" = "#cc7818",
  "Haematoendothelial_progenitors" = "#FBBE92",
  "Endothelium" = "#ff891c",
  "Blood_progenitors" = "#c9a997",
  "Blood_progenitors_1" = "#f9decf",
  "Blood_progenitors_2" = "#c9a997",
  "Erythroid" = "#EF4E22",
  "Erythroid1" = "#C72228",
  "Erythroid2" = "#f79083",
  "Erythroid3" = "#EF4E22",
  "NMP" = "#8EC792",
  "Rostral_neurectoderm" = "#65A83E",
  "Caudal_neurectoderm" = "#354E23",
  "Neural_crest" = "#C3C388",
  "Forebrain_Midbrain_Hindbrain" = "#647a4f",
  "Spinal_cord" = "#CDE088",
  "Surface_ectoderm" = "#f7f79e",
  "Visceral_endoderm" = "#F6BFCB",
  "ExE_endoderm" = "#7F6874",
  "ExE_ectoderm" = "#989898",
  "Parietal_endoderm" = "#1A1A1A",
  
  # Additional
  "Erythroid" = "#EF4E22",
  "Blood_progenitors" = "#c9a997"
)

opts$samples <- c(
  # "E125_DNMT3A_HET_A_L001",
  # "E125_DNMT3A_HET_A_L003",
  # "E125_DNMT3A_KO_B_L002",
  # "E125_DNMT3A_KO_E_L004",
  # "A_E12_5_D3a_Het_L001",
  # "B_E12_5_D3a_KO_L002",
  "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001",
  "SIGAB6_E85_3_Dnmt3aWT_Dnmt3b_WT_L002",
  "SIGAC6_E85_5_Dnmt3aKO_Dnmt3b_Het_L003",
  "SIGAD6_E85_8_Dnmt3aHet_Dnmt3b_KO_L004",
  "15_E8_5_D3A_WT_D3B_WT_L007",
  "17_E8_5_D3A_KO_D3B_WT_L008",
  "2_E8_5_D3A_WT_D3B_KO_L003",
  "3_E8_5_D3A_HET_D3B_WT_L004",
  "7_E8_5_D3A_WT_D3B_KO_L005",
  "8_E8_5_D3A_KO_D3B_KO_L006",
  "E8_5_Dnmt1_KO_male_SIGAC8_L001",
  "E8_5_Dnmt1_KO_male_SIGAD8_L002",
  "E8_5_Dnmt1_KO_male_SIGAE8_L003",
  "E8_5_Dnmt1_WT_female_SIGAB8_L004",
  "E8_5_Dnmt1_WT_female_SIGAF8_L005",
  "E8_5_Dnmt3ab_WT_female_SIGAA8_L006",
  "SIGAH10_Dnmt3ab_WT_L002",
  "SIGAH11_Dnmt3ab_WT_L003",
  "SIGAH9_Dnmt3a_KO_Dnmt3b_Het_L001",
  "SIGAG5_9_dnmt3ab_DKO_L005"
)

opts$sample2class <- c(
  # "E125_DNMT3A_HET_A_L001" = "E12.5_Dnmt3aWT_Dnmt3bHET",
  # "E125_DNMT3A_HET_A_L003" = "E12.5_Dnmt3aWT_Dnmt3bHET",
  # "E125_DNMT3A_KO_B_L002" = "E12.5_Dnmt3aWT_Dnmt3bKO",
  # "E125_DNMT3A_KO_E_L004" = "E12.5_Dnmt3aWT_Dnmt3bKO",
  # "A_E12_5_D3a_Het_L001" = "E12.5_Dnmt3aHET_Dnmt3bWT",
  # "B_E12_5_D3a_KO_L002"  = "E12.5_Dnmt3aKO_Dnmt3bWT ",
  "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001" = "E8.5_Dnmt3aKO_Dnmt3bWT",
  "SIGAB6_E85_3_Dnmt3aWT_Dnmt3b_WT_L002" = "E8.5_WT",
  "SIGAC6_E85_5_Dnmt3aKO_Dnmt3b_Het_L003" = "E8.5_Dnmt3aKO_Dnmt3bHET",
  "SIGAD6_E85_8_Dnmt3aHet_Dnmt3b_KO_L004" = "E8.5_Dnmt3aHET_Dnmt3bKO",
  "15_E8_5_D3A_WT_D3B_WT_L007" = "E8.5_WT",
  "17_E8_5_D3A_KO_D3B_WT_L008" = "E8.5_Dnmt3aKO_Dnmt3bWT",
  "2_E8_5_D3A_WT_D3B_KO_L003" = "E8.5_Dnmt3aWT_Dnmt3bKO",
  "3_E8_5_D3A_HET_D3B_WT_L004" = "E8.5_Dnmt3aHET_Dnmt3bWT",
  "7_E8_5_D3A_WT_D3B_KO_L005" = "E8.5_Dnmt3aWT_Dnmt3bKO",
  "8_E8_5_D3A_KO_D3B_KO_L006" = "E8.5_Dnmt3aKO_Dnmt3bKO",
  "E8_5_Dnmt1_KO_male_SIGAC8_L001" = "E8.5_Dnmt1KO",
  "E8_5_Dnmt1_KO_male_SIGAD8_L002" = "E8.5_Dnmt1KO",
  "E8_5_Dnmt1_KO_male_SIGAE8_L003" = "E8.5_Dnmt1KO",
  "E8_5_Dnmt1_WT_female_SIGAB8_L004" = "E8.5_WT",
  "E8_5_Dnmt1_WT_female_SIGAF8_L005" = "E8.5_WT",
  "E8_5_Dnmt3ab_WT_female_SIGAA8_L006" = "E8.5_WT",
  "SIGAH10_Dnmt3ab_WT_L002" = "E8.5_WT",
  "SIGAH11_Dnmt3ab_WT_L003" = "E8.5_WT",
  "SIGAH9_Dnmt3a_KO_Dnmt3b_Het_L001" = "E8.5_Dnmt3aKO_Dnmt3bHET",
  "SIGAG5_9_dnmt3ab_DKO_L005" = "E8.5_Dnmt3aKO_Dnmt3bKO"
)

opts$stages <- c(
  # "E12.5",
  "E8.5"
)

opts$classes <- c(
  # "E12.5_Dnmt3aWT_Dnmt3bHET",
  # "E12.5_Dnmt3aWT_Dnmt3bKO",
  # "E12.5_Dnmt3aHET_Dnmt3bWT",
  # "E12.5_Dnmt3aKO_Dnmt3bWT",
  "E8.5_Dnmt3aKO_Dnmt3bWT", 
  "E8.5_WT", 
  "E8.5_Dnmt3aHET_Dnmt3bKO", 
  "E8.5_Dnmt3aHET_Dnmt3bWT", 
  "E8.5_Dnmt3aKO_Dnmt3bHET", 
  "E8.5_Dnmt3aKO_Dnmt3bKO", 
  "E8.5_Dnmt3aWT_Dnmt3bKO",
  "E8.5_Dnmt1KO"
)

opts$sample2alias <- c(
  # "E125_DNMT3A_HET_A_L001" = "E12.5_Dnmt3aWT_Dnmt3bHET_1",
  # "E125_DNMT3A_HET_A_L003" = "E12.5_Dnmt3aWT_Dnmt3bHET_2",
  # "E125_DNMT3A_KO_B_L002" = "E12.5_Dnmt3aWT_Dnmt3bKO_1",
  # "E125_DNMT3A_KO_E_L004" = "E12.5_Dnmt3aWT_Dnmt3bKO_2",
  # "A_E12_5_D3a_Het_L001" = "E12.5_Dnmt3aHET_Dnmt3bWT",
  # "B_E12_5_D3a_KO_L002"  = "E12.5_Dnmt3aKO_Dnmt3bWT ",
  "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001" = "E8.5_Dnmt3aKO_Dnmt3bWT_1",
  "17_E8_5_D3A_KO_D3B_WT_L008" = "E8.5_Dnmt3aKO_Dnmt3bWT_2",
  "SIGAB6_E85_3_Dnmt3aWT_Dnmt3b_WT_L002" = "E8.5_WT_1",
  "15_E8_5_D3A_WT_D3B_WT_L007" = "E8.5_WT_2",
  "E8_5_Dnmt1_WT_female_SIGAB8_L004" = "E8.5_WT_3",
  "E8_5_Dnmt1_WT_female_SIGAF8_L005" = "E8.5_WT_4",
  "E8_5_Dnmt3ab_WT_female_SIGAA8_L006" = "E8.5_WT_5",
  "SIGAH10_Dnmt3ab_WT_L002" = "E8.5_WT_6",
  "SIGAH11_Dnmt3ab_WT_L003" = "E8.5_WT_7",
  "SIGAC6_E85_5_Dnmt3aKO_Dnmt3b_Het_L003" = "E8.5_Dnmt3aKO_Dnmt3bHET_1",
  "SIGAH9_Dnmt3a_KO_Dnmt3b_Het_L001" = "E8.5_Dnmt3aKO_Dnmt3bHET_2",
  "SIGAD6_E85_8_Dnmt3aHet_Dnmt3b_KO_L004" = "E8.5_Dnmt3aHET_Dnmt3bKO",
  "2_E8_5_D3A_WT_D3B_KO_L003" = "E8.5_Dnmt3aWT_Dnmt3bKO_1",
  "7_E8_5_D3A_WT_D3B_KO_L005" = "E8.5_Dnmt3aWT_Dnmt3bKO_2",
  "3_E8_5_D3A_HET_D3B_WT_L004" = "E8.5_Dnmt3aHET_Dnmt3bWT",
  "8_E8_5_D3A_KO_D3B_KO_L006" = "E8.5_Dnmt3aKO_Dnmt3bKO_1",
  "SIGAG5_9_dnmt3ab_DKO_L005" = "E8.5_Dnmt3aKO_Dnmt3bKO_2",
  "E8_5_Dnmt1_KO_male_SIGAC8_L001" = "E8.5_Dnmt1KO_1",
  "E8_5_Dnmt1_KO_male_SIGAD8_L002" = "E8.5_Dnmt1KO_2",
  "E8_5_Dnmt1_KO_male_SIGAE8_L003" = "E8.5_Dnmt1KO_3"
)


opts$aliases <- c(
  # "E12.5_Dnmt3aHET_Dnmt3bWT",
  # "E12.5_Dnmt3aKO_Dnmt3bWT",
  # "E12.5_Dnmt3aWT_Dnmt3bHET_1",
  # "E12.5_Dnmt3aWT_Dnmt3bHET_2",
  # "E12.5_Dnmt3aWT_Dnmt3bKO_1",
  # "E12.5_Dnmt3aWT_Dnmt3bKO_2",
  "E8.5_Dnmt1KO_1",
  "E8.5_Dnmt1KO_2",
  "E8.5_Dnmt1KO_3",
  "E8.5_Dnmt3aHET_Dnmt3bKO",
  "E8.5_Dnmt3aHET_Dnmt3bWT",
  "E8.5_Dnmt3aKO_Dnmt3bHET_1",
  "E8.5_Dnmt3aKO_Dnmt3bHET_2",
  "E8.5_Dnmt3aKO_Dnmt3bKO_1",
  "E8.5_Dnmt3aKO_Dnmt3bKO_2",
  "E8.5_Dnmt3aKO_Dnmt3bWT_1",
  "E8.5_Dnmt3aKO_Dnmt3bWT_2",
  "E8.5_Dnmt3aWT_Dnmt3bKO_1",
  "E8.5_Dnmt3aWT_Dnmt3bKO_2",
  "E8.5_WT_1",
  "E8.5_WT_2",
  "E8.5_WT_3",
  "E8.5_WT_4",
  "E8.5_WT_5",
  "E8.5_WT_6",
  "E8.5_WT_7"
)

opts$ExE.celltypes <- c(
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm",
  "Parietal_endoderm"
)

##########################
## Load sample metadata ##
##########################

# sample_metadata <- fread(io$metadata) %>% .[pass_QC==T]# %>%
  # .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped," ","_")] %>%
  # .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,"/","_")]

  