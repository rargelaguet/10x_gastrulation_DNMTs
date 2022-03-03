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
io$seurat <- paste0(io$basedir,"/processed/seurat.rds")
io$sce <- paste0(io$basedir,"/processed/SingleCellExperiment.rds")
# io$sex <- paste0(io$basedir,"/results/sex/sex_assignment.txt.gz")

# Atlas information
io$atlas.metadata <- paste0(io$atlas.basedir,"/sample_metadata.txt.gz")
io$atlas.marker_genes <- paste0(io$atlas.basedir,"/results/marker_genes/all_stages/marker_genes.txt.gz")
io$atlas.differential <- paste0(io$atlas.basedir,"/results/differential")
# io$atlas.average_expression_per_celltype <- paste0(io$atlas.basedir,"/results/marker_genes/avg_expr_per_celltype_and_gene.txt.gz")
io$atlas.sce <- paste0(io$atlas.basedir,"/processed/SingleCellExperiment.rds")

# paga
io$paga.connectivity <- file.path(io$atlas.basedir,"results/paga/paga_connectivity.csv")
io$paga.coordinates <- file.path(io$atlas.basedir,"results/paga/paga_coordinates.csv")


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
  # "3_E8_5_D3A_HET_D3B_WT_L004",
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
  "SIGAG5_9_dnmt3ab_DKO_L005",
  "Dnmt1_E8.5_embryo10_Grosswendt2020",
  "Dnmt1_E8.5_embryo11_Grosswendt2020",
  "Dnmt1_E8.5_embryo12_Grosswendt2020",
  "Dnmt1_E8.5_embryo1_Grosswendt2020",
  "Dnmt1_E8.5_embryo2_Grosswendt2020",
  "Dnmt1_E8.5_embryo3_Grosswendt2020",
  "Dnmt1_E8.5_embryo4_Grosswendt2020",
  "Dnmt1_E8.5_embryo5_Grosswendt2020",
  "Dnmt1_E8.5_embryo6_Grosswendt2020",
  "Dnmt1_E8.5_embryo7_Grosswendt2020",
  "Dnmt1_E8.5_embryo8_Grosswendt2020",
  "Dnmt1_E8.5_embryo9_Grosswendt2020",
  "Dnmt3a_E8.5_embryo10_Grosswendt2020",
  "Dnmt3a_E8.5_embryo1_Grosswendt2020",
  "Dnmt3a_E8.5_embryo2_Grosswendt2020",
  "Dnmt3a_E8.5_embryo3_Grosswendt2020",
  "Dnmt3a_E8.5_embryo4_Grosswendt2020",
  "Dnmt3a_E8.5_embryo5_Grosswendt2020",
  "Dnmt3a_E8.5_embryo6_Grosswendt2020",
  "Dnmt3a_E8.5_embryo7_Grosswendt2020",
  "Dnmt3a_E8.5_embryo8_Grosswendt2020",
  "Dnmt3a_E8.5_embryo9_Grosswendt2020",
  "Dnmt3b_E8.5_embryo1_Grosswendt2020",
  "Dnmt3b_E8.5_embryo2_Grosswendt2020",
  "Dnmt3b_E8.5_embryo3_Grosswendt2020",
  "Dnmt3b_E8.5_embryo4_Grosswendt2020",
  "Dnmt3b_E8.5_embryo5_Grosswendt2020",
  "Dnmt3b_E8.5_embryo6_Grosswendt2020",
  "Dnmt3b_E8.5_embryo7_Grosswendt2020", 
  "Dnmt3b_E8.5_embryo8_Grosswendt2020",
  "WT_E8.5_embryo1_Grosswendt2020",
  "WT_E8.5_embryo2_Grosswendt2020",
  "WT_E8.5_embryo3_Grosswendt2020",
  "WT_E8.5_embryo4_Grosswendt2020",
  "WT_E8.5_embryo5_Grosswendt2020",
  "WT_E8.5_embryo6_Grosswendt2020",
  "WT_E8.5_embryo7_Grosswendt2020",
  "WT_E8.5_embryo8_Grosswendt2020",
  "WT_E8.5_embryo9_Grosswendt2020",
  "WT_E8.5_embryo10_Grosswendt2020"
)

opts$sample2class <- c(
  # "E125_DNMT3A_HET_A_L001" = "E12.5_Dnmt3aWT_Dnmt3bHET",
  # "E125_DNMT3A_HET_A_L003" = "E12.5_Dnmt3aWT_Dnmt3bHET",
  # "E125_DNMT3A_KO_B_L002" = "E12.5_Dnmt3b_KO",
  # "E125_DNMT3A_KO_E_L004" = "E12.5_Dnmt3b_KO",
  # "A_E12_5_D3a_Het_L001" = "E12.5_Dnmt3a_HET_Dnmt3b_WT",
  # "B_E12_5_D3a_KO_L002"  = "E12.5_Dnmt3a_KO ",
  "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001" = "Dnmt3a_KO",
  "SIGAB6_E85_3_Dnmt3aWT_Dnmt3b_WT_L002" = "WT",
  # "SIGAC6_E85_5_Dnmt3aKO_Dnmt3b_Het_L003" = "Dnmt3a_KO_Dnmt3b_HET",
  "SIGAC6_E85_5_Dnmt3aKO_Dnmt3b_Het_L003" = "Dnmt3a_KO",
  # "SIGAD6_E85_8_Dnmt3aHet_Dnmt3b_KO_L004" = "Dnmt3a_HET_Dnmt3b_KO",
  "SIGAD6_E85_8_Dnmt3aHet_Dnmt3b_KO_L004" = "Dnmt3b_KO",
  "15_E8_5_D3A_WT_D3B_WT_L007" = "WT",
  "17_E8_5_D3A_KO_D3B_WT_L008" = "Dnmt3a_KO",
  "2_E8_5_D3A_WT_D3B_KO_L003" = "Dnmt3b_KO",
  # "3_E8_5_D3A_HET_D3B_WT_L004" = "Dnmt3a_HET_Dnmt3b_WT",
  "7_E8_5_D3A_WT_D3B_KO_L005" = "Dnmt3b_KO",
  "8_E8_5_D3A_KO_D3B_KO_L006" = "Dnmt3ab_KO",
  "E8_5_Dnmt1_KO_male_SIGAC8_L001" = "Dnmt1_KO",
  "E8_5_Dnmt1_KO_male_SIGAD8_L002" = "Dnmt1_KO",
  "E8_5_Dnmt1_KO_male_SIGAE8_L003" = "Dnmt1_KO",
  "E8_5_Dnmt1_WT_female_SIGAB8_L004" = "WT",
  "E8_5_Dnmt1_WT_female_SIGAF8_L005" = "WT",
  "E8_5_Dnmt3ab_WT_female_SIGAA8_L006" = "WT",
  "SIGAH10_Dnmt3ab_WT_L002" = "WT",
  "SIGAH11_Dnmt3ab_WT_L003" = "WT",
  # "SIGAH9_Dnmt3a_KO_Dnmt3b_Het_L001" = "Dnmt3a_KO_Dnmt3b_HET",
  "SIGAH9_Dnmt3a_KO_Dnmt3b_Het_L001" = "Dnmt3a_KO",
  # "SIGAG5_9_dnmt3ab_DKO_L005" = "Dnmt3ab_KO",
  "SIGAG5_9_dnmt3ab_DKO_L005" = "Dnmt3b_KO",
  "Dnmt1_E8.5_embryo10_Grosswendt2020" = "Dnmt1_KO",
  "Dnmt1_E8.5_embryo11_Grosswendt2020" = "Dnmt1_KO",
  "Dnmt1_E8.5_embryo12_Grosswendt2020" = "Dnmt1_KO",
  "Dnmt1_E8.5_embryo1_Grosswendt2020" = "Dnmt1_KO",
  "Dnmt1_E8.5_embryo2_Grosswendt2020" = "Dnmt1_KO",
  "Dnmt1_E8.5_embryo3_Grosswendt2020" = "Dnmt1_KO",
  "Dnmt1_E8.5_embryo4_Grosswendt2020" = "Dnmt1_KO",
  "Dnmt1_E8.5_embryo5_Grosswendt2020" = "Dnmt1_KO",
  "Dnmt1_E8.5_embryo6_Grosswendt2020" = "Dnmt1_KO",
  "Dnmt1_E8.5_embryo7_Grosswendt2020" = "Dnmt1_KO",
  "Dnmt1_E8.5_embryo8_Grosswendt2020" = "Dnmt1_KO",
  "Dnmt1_E8.5_embryo9_Grosswendt2020" = "Dnmt1_KO",
  "Dnmt3a_E8.5_embryo10_Grosswendt2020" = "Dnmt3a_KO",
  "Dnmt3a_E8.5_embryo1_Grosswendt2020" = "Dnmt3a_KO",
  "Dnmt3a_E8.5_embryo2_Grosswendt2020" = "Dnmt3a_KO",
  "Dnmt3a_E8.5_embryo3_Grosswendt2020" = "Dnmt3a_KO",
  "Dnmt3a_E8.5_embryo4_Grosswendt2020" = "Dnmt3a_KO",
  "Dnmt3a_E8.5_embryo5_Grosswendt2020" = "Dnmt3a_KO",
  "Dnmt3a_E8.5_embryo6_Grosswendt2020" = "Dnmt3a_KO",
  "Dnmt3a_E8.5_embryo7_Grosswendt2020" = "Dnmt3a_KO",
  "Dnmt3a_E8.5_embryo8_Grosswendt2020" = "Dnmt3a_KO",
  "Dnmt3a_E8.5_embryo9_Grosswendt2020" = "Dnmt3a_KO",
  "Dnmt3b_E8.5_embryo1_Grosswendt2020" = "Dnmt3b_KO",
  "Dnmt3b_E8.5_embryo2_Grosswendt2020" = "Dnmt3b_KO",
  "Dnmt3b_E8.5_embryo3_Grosswendt2020" = "Dnmt3b_KO",
  "Dnmt3b_E8.5_embryo4_Grosswendt2020" = "Dnmt3b_KO",
  "Dnmt3b_E8.5_embryo5_Grosswendt2020" = "Dnmt3b_KO",
  "Dnmt3b_E8.5_embryo6_Grosswendt2020" = "Dnmt3b_KO",
  "Dnmt3b_E8.5_embryo7_Grosswendt2020" = "Dnmt3b_KO",
  "Dnmt3b_E8.5_embryo8_Grosswendt2020" = "Dnmt3b_KO",
  "WT_E8.5_embryo1_Grosswendt2020" = "WT",
  "WT_E8.5_embryo2_Grosswendt2020" = "WT",
  "WT_E8.5_embryo3_Grosswendt2020" = "WT",
  "WT_E8.5_embryo4_Grosswendt2020" = "WT",
  "WT_E8.5_embryo5_Grosswendt2020" = "WT",
  "WT_E8.5_embryo6_Grosswendt2020" = "WT",
  "WT_E8.5_embryo7_Grosswendt2020" = "WT",
  "WT_E8.5_embryo8_Grosswendt2020" = "WT",
  "WT_E8.5_embryo9_Grosswendt2020" = "WT",
  "WT_E8.5_embryo10_Grosswendt2020" = "WT"
)

opts$stages <- c(
  # "E12.5",
  "E8.5"
)

opts$atlas.stages <- c(
  "E6.5",
  "E6.75",
  "E7.0",
  "E7.25",
  "E7.5",
  "E7.75",
  "E8.0",
  "E8.25",
  "E8.5"
  # "mixed_gastrulation"
)

opts$classes <- c(
  # "E12.5_Dnmt3aWT_Dnmt3bHET",
  # "E12.5_Dnmt3b_KO",
  # "E12.5_Dnmt3a_HET_Dnmt3b_WT",
  # "E12.5_Dnmt3a_KO",
  "WT", 
  "Dnmt3a_KO", 
  # "Dnmt3a_HET_Dnmt3b_KO", 
  # "Dnmt3a_HET_Dnmt3b_WT", 
  # "Dnmt3a_KO_Dnmt3b_HET", 
  "Dnmt3b_KO",
  "Dnmt1_KO",
  "Dnmt3ab_KO"
)


opts$classes.colors <- c(
  "WT" = "#ffffb3", 
  "Dnmt3a_KO" = "#8dd3c7", 
  "Dnmt3b_KO" = "#fb8072", 
  "Dnmt1_KO" = "#80b1d3",
  "Dnmt3ab_KO" = "#bebada"
)

opts$stage.colors = c(
  "E8.5" = "#440154FF",
  "E8.25" = "#472D7BFF",
  "E8.0" = "#3B528BFF",
  "E7.75" = "#2C728EFF",
  "E7.5" = "#21908CFF",
  "E7.25" = "#27AD81FF",
  "E7.0" = "#5DC863FF",
  "E6.75" = "#AADC32FF",
  "E6.5" = "#FDE725FF"
)

opts$sample2alias <- c(
  # "E125_DNMT3A_HET_A_L001" = "E12.5_Dnmt3aWT_Dnmt3bHET_1",
  # "E125_DNMT3A_HET_A_L003" = "E12.5_Dnmt3aWT_Dnmt3bHET_2",
  # "E125_DNMT3A_KO_B_L002" = "E12.5_Dnmt3b_KO_1",
  # "E125_DNMT3A_KO_E_L004" = "E12.5_Dnmt3b_KO_2",
  # "A_E12_5_D3a_Het_L001" = "E12.5_Dnmt3a_HET_Dnmt3b_WT",
  # "B_E12_5_D3a_KO_L002"  = "E12.5_Dnmt3a_KO ",
  "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001" = "Dnmt3a_KO_1",
  "17_E8_5_D3A_KO_D3B_WT_L008" = "Dnmt3a_KO_2",
  "SIGAB6_E85_3_Dnmt3aWT_Dnmt3b_WT_L002" = "WT_1",
  "15_E8_5_D3A_WT_D3B_WT_L007" = "WT_2",
  "E8_5_Dnmt1_WT_female_SIGAB8_L004" = "WT_3",
  "E8_5_Dnmt1_WT_female_SIGAF8_L005" = "WT_4",
  "E8_5_Dnmt3ab_WT_female_SIGAA8_L006" = "WT_5",
  "SIGAH10_Dnmt3ab_WT_L002" = "WT_6",
  "SIGAH11_Dnmt3ab_WT_L003" = "WT_7",
  # "SIGAC6_E85_5_Dnmt3aKO_Dnmt3b_Het_L003" = "Dnmt3a_KO_Dnmt3b_HET_1",
  "SIGAC6_E85_5_Dnmt3aKO_Dnmt3b_Het_L003" = "Dnmt3a_KO_13",
  # "SIGAH9_Dnmt3a_KO_Dnmt3b_Het_L001" = "Dnmt3a_KO_Dnmt3b_HET_2",
  "SIGAH9_Dnmt3a_KO_Dnmt3b_Het_L001" = "Dnmt3a_KO_14",
  # "SIGAD6_E85_8_Dnmt3aHet_Dnmt3b_KO_L004" = "Dnmt3a_HET_Dnmt3b_KO",
  "SIGAD6_E85_8_Dnmt3aHet_Dnmt3b_KO_L004" = "Dnmt3b_KO_11",
  "2_E8_5_D3A_WT_D3B_KO_L003" = "Dnmt3b_KO_1",
  "7_E8_5_D3A_WT_D3B_KO_L005" = "Dnmt3b_KO_2",
  # "3_E8_5_D3A_HET_D3B_WT_L004" = "Dnmt3a_HET_Dnmt3b_WT",
  "8_E8_5_D3A_KO_D3B_KO_L006" = "Dnmt3ab_KO_1",
  # "SIGAG5_9_dnmt3ab_DKO_L005" = "Dnmt3ab_KO_2",
  "SIGAG5_9_dnmt3ab_DKO_L005" = "Dnmt3b_KO_12",
  "E8_5_Dnmt1_KO_male_SIGAC8_L001" = "Dnmt1_KO_1",
  "E8_5_Dnmt1_KO_male_SIGAD8_L002" = "Dnmt1_KO_2",
  "E8_5_Dnmt1_KO_male_SIGAE8_L003" = "Dnmt1_KO_3",
  "Dnmt1_E8.5_embryo1_Grosswendt2020" = "Dnmt1_KO_4",
  "Dnmt1_E8.5_embryo2_Grosswendt2020" = "Dnmt1_KO_5",
  "Dnmt1_E8.5_embryo3_Grosswendt2020" = "Dnmt1_KO_6",
  "Dnmt1_E8.5_embryo4_Grosswendt2020" = "Dnmt1_KO_7",
  "Dnmt1_E8.5_embryo5_Grosswendt2020" = "Dnmt1_KO_8",
  "Dnmt1_E8.5_embryo6_Grosswendt2020" = "Dnmt1_KO_9",
  "Dnmt1_E8.5_embryo7_Grosswendt2020" = "Dnmt1_KO_10",
  "Dnmt1_E8.5_embryo8_Grosswendt2020" = "Dnmt1_KO_11",
  "Dnmt1_E8.5_embryo9_Grosswendt2020" = "Dnmt1_KO_12",
  "Dnmt1_E8.5_embryo10_Grosswendt2020" = "Dnmt1_KO_13",
  "Dnmt1_E8.5_embryo11_Grosswendt2020" = "Dnmt1_KO_14",
  "Dnmt1_E8.5_embryo12_Grosswendt2020" = "Dnmt1_KO_15",
  "Dnmt3a_E8.5_embryo10_Grosswendt2020" = "Dnmt3a_KO_3",
  "Dnmt3a_E8.5_embryo1_Grosswendt2020" = "Dnmt3a_KO_4",
  "Dnmt3a_E8.5_embryo2_Grosswendt2020" = "Dnmt3a_KO_5",
  "Dnmt3a_E8.5_embryo3_Grosswendt2020" = "Dnmt3a_KO_6",
  "Dnmt3a_E8.5_embryo4_Grosswendt2020" = "Dnmt3a_KO_7",
  "Dnmt3a_E8.5_embryo5_Grosswendt2020" = "Dnmt3a_KO_8",
  "Dnmt3a_E8.5_embryo6_Grosswendt2020" = "Dnmt3a_KO_9",
  "Dnmt3a_E8.5_embryo7_Grosswendt2020" = "Dnmt3a_KO_10",
  "Dnmt3a_E8.5_embryo8_Grosswendt2020" = "Dnmt3a_KO_11",
  "Dnmt3a_E8.5_embryo9_Grosswendt2020" = "Dnmt3a_KO_12",
  "Dnmt3b_E8.5_embryo1_Grosswendt2020" = "Dnmt3b_KO_3",
  "Dnmt3b_E8.5_embryo2_Grosswendt2020" = "Dnmt3b_KO_4",
  "Dnmt3b_E8.5_embryo3_Grosswendt2020" = "Dnmt3b_KO_5",
  "Dnmt3b_E8.5_embryo4_Grosswendt2020" = "Dnmt3b_KO_6",
  "Dnmt3b_E8.5_embryo5_Grosswendt2020" = "Dnmt3b_KO_7",
  "Dnmt3b_E8.5_embryo6_Grosswendt2020" = "Dnmt3b_KO_8",
  "Dnmt3b_E8.5_embryo7_Grosswendt2020" = "Dnmt3b_KO_9",
  "Dnmt3b_E8.5_embryo8_Grosswendt2020" = "Dnmt3b_KO_10",
  "WT_E8.5_embryo1_Grosswendt2020" = "WT_8",
  "WT_E8.5_embryo2_Grosswendt2020" = "WT_9",
  "WT_E8.5_embryo3_Grosswendt2020" = "WT_10",
  "WT_E8.5_embryo4_Grosswendt2020" = "WT_11",
  "WT_E8.5_embryo5_Grosswendt2020" = "WT_12",
  "WT_E8.5_embryo6_Grosswendt2020" = "WT_13",
  "WT_E8.5_embryo7_Grosswendt2020" = "WT_14",
  "WT_E8.5_embryo8_Grosswendt2020" = "WT_15",
  "WT_E8.5_embryo9_Grosswendt2020" = "WT_16",
  "WT_E8.5_embryo10_Grosswendt2020" = "WT_17"
)

# opts$aliases <- c(
#   # "E12.5_Dnmt3a_HET_Dnmt3b_WT",
#   # "E12.5_Dnmt3a_KO",
#   # "E12.5_Dnmt3aWT_Dnmt3bHET_1",
#   # "E12.5_Dnmt3aWT_Dnmt3bHET_2",
#   # "E12.5_Dnmt3b_KO_1",
#   # "E12.5_Dnmt3b_KO_2",
#   "Dnmt1_KO_1",
#   "Dnmt1_KO_2",
#   "Dnmt1_KO_3",
#   "Dnmt3a_HET_Dnmt3b_KO",
#   "Dnmt3a_HET_Dnmt3b_WT",
#   "Dnmt3a_KO_Dnmt3b_HET_1",
#   "Dnmt3a_KO_Dnmt3b_HET_2",
#   "Dnmt3ab_KO_1",
#   "Dnmt3ab_KO_2",
#   "Dnmt3a_KO_1",
#   "Dnmt3a_KO_2",
#   "Dnmt3b_KO_1",
#   "Dnmt3b_KO_2",
#   "WT_1",
#   "WT_2",
#   "WT_3",
#   "WT_4",
#   "WT_5",
#   "WT_6",
#   "WT_7",
# )

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

  