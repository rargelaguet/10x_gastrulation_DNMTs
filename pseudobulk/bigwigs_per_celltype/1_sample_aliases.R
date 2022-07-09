source(here::here("settings.R"))

here::i_am("pseudobulk/biwigs_per_celltype/1_sample_aliases.R")


# this script loads our metadata file then generates an alias file to match barcodes to e.g. celltypes
# this is then fed into the subsequent samtools scripts to split and merge 10x bam files by e.g. celltypes
 
######################
## Define arguments ##
######################

p <- ArgumentParser(description='')

p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--outfile',    type="character",    help='Alias file to output')
p$add_argument('--group_by',    type="character",    help='Metadata column to group cells by')
p$add_argument('--indir',      type="character",    help='10X BAM file parent directory')
p$add_argument('--outdir',      type="character",    help='Output directory')

args <- p$parse_args(commandArgs(TRUE))



#####################
## Define settings ##
#####################

# match sample names to 10x output folders:

folders <- c(
  "CellRanger6_L001_SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT" = "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001" , 
  "CellRanger6_L006_8_E8_5_D3A_KO_D3B_KO" = "8_E8_5_D3A_KO_D3B_KO_L006",  
  "CellRanger6_L001_SIGAH9_Dnmt3a_KO_Dnmt3b_Het"   = "SIGAH9_Dnmt3a_KO_Dnmt3b_Het_L001",    
  "CellRanger6_L007_15_E8_5_D3A_WT_D3B_WT" = "15_E8_5_D3A_WT_D3B_WT_L007",
  "CellRanger6_L002_SIGAB6_E85_3_Dnmt3aWT_Dnmt3b_WT"   = "SIGAB6_E85_3_Dnmt3aWT_Dnmt3b_WT_L002",
  "CellRanger6_L008_17_E8_5_D3A_KO_D3B_WT" = "17_E8_5_D3A_KO_D3B_WT_L008",
  "CellRanger6_L002_SIGAH10_Dnmt3ab_WT"       = "SIGAH10_Dnmt3ab_WT_L002",         
  "CellRanger6_SIGAG5_9_dnmt3ab_DKO" = "SIGAG5_9_dnmt3ab_DKO_L005",
  "CellRanger6_L003_2_E8_5_D3A_WT_D3B_KO"   = "2_E8_5_D3A_WT_D3B_KO_L003",
  "E8_5_Dnmt1_KO_male_SIGAC8" = "E8_5_Dnmt1_KO_male_SIGAC8_L001",
  "CellRanger6_L003_SIGAC6_E85_5_Dnmt3aKO_Dnmt3b_Het"  = "SIGAC6_E85_5_Dnmt3aKO_Dnmt3b_Het_L003",
  "E8_5_Dnmt1_KO_male_SIGAD8" = "E8_5_Dnmt1_KO_male_SIGAD8_L002",
  "CellRanger6_L003_SIGAH11_Dnmt3ab_WT" = "SIGAH11_Dnmt3ab_WT_L003",
  "E8_5_Dnmt1_KO_male_SIGAE8" = "E8_5_Dnmt1_KO_male_SIGAE8_L003",
  #CellRanger6_L004_3_E8_5_D3A_HET_D3B_WT             
  "E8_5_Dnmt1_WT_female_SIGAB8" = "E8_5_Dnmt1_WT_female_SIGAB8_L004",
  "CellRanger6_L004_SIGAD6_E85_8_Dnmt3aHet_Dnmt3b_KO"  = "SIGAD6_E85_8_Dnmt3aHet_Dnmt3b_KO_L004",
  "E8_5_Dnmt1_WT_female_SIGAF8" = "E8_5_Dnmt1_WT_female_SIGAF8_L005",
  "CellRanger6_L005_7_E8_5_D3A_WT_D3B_KO"   =  "7_E8_5_D3A_WT_D3B_KO_L005" , 
  "E8_5_Dnmt3ab_WT_female_SIGAA8" = "E8_5_Dnmt3ab_WT_female_SIGAA8_L006" 
)


# ## START TEST ##
# args$metadata <- file.path(io$basedir,"/sample_metadata_after_mapping.txt.gz")
# 
# samples <- fread(args$metadata) %>% 
#   .[dataset=="KO", unique(sample)]
# length(samples)
# 
# args$group_by <- "class_sample_celltype"
# args$sample <- samples[12]

# #args$outdir <- file.path(io$basedir,"results/pseudobulk/bam")
# args$outdir <- "/bi/scratch/Stephen_Clark/dnmt_10x_bam_files_for_genotyping/pseudobulk"
# args$indir  <- "/bi/scratch/Stephen_Clark/dnmt_10x_bam_files_for_genotyping"
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings=F)

# Options
# opts$rename_celltypes <- c(
#   "Erythroid3" = "Erythroid",
#   "Erythroid2" = "Erythroid",
#   "Erythroid1" = "Erythroid",
#   "Blood_progenitors_1" = "Blood_progenitors",
#   "Blood_progenitors_2" = "Blood_progenitors",
#   "Anterior_Primitive_Streak" = "Primitive_Streak"
# )


###############
## Load data ##
###############

# Load cell metadata
sample_metadata <- fread(args$metadata) %>%
  .[,dataset:=ifelse(grepl("Grosswendt",sample),"CRISPR","KO")] %>%
  # .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,opts$rename_celltypes)] %>%
  .[,class_celltype:=sprintf("%s_%s",class,celltype.mapped)] %>%
  .[,class_celltype_dataset:=sprintf("%s_%s_%s",class,celltype.mapped,dataset)] %>%
  .[,class_sample_celltype:=sprintf("%s_%s_%s",class,sample,celltype.mapped)] %>%
  # .[,class_celltype_dataset:=sprintf("%s-%s-%s",class,celltype.mapped,dataset)] %>%
  .[pass_rnaQC==TRUE & !is.na(eval(as.name(args$group_by)))]

table(sample_metadata[[args$group_by]])
samples <- sample_metadata[, unique(sample)]

# fitler to exclude crispr data
sample_metadata <- sample_metadata[dataset == "KO"]

# format barcode so it matches the BAM 
if (sample_metadata[1, !barcode %like% "CB:Z:"]){
  sample_metadata[, barcode := paste0("CB:Z:", barcode)]
}


# add in bam file location

sample_metadata <- merge(
  sample_metadata, 
  data.table(sample = folders, folder = names(folders)),
  by = "sample"
)

sample_metadata[, path := paste0(args$indir, "/", folder, "/possorted_genome_bam.bam")]


aliases <- sample_metadata[, .(
  sample, 
  bam_path = path,
  barcode,
  split_to = paste0(args$outdir, "/", sample, "/", sample, "_", get(args$group_by), ".bam"),
  merge_to = paste0(args$outdir, "/", get(args$group_by), ".bam"),
  outdir = args$outdir
)]


fwrite(aliases, args$outfile, sep = "\t", na = "NA", quote = FALSE)




