
################
## Define I/O ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
  io$script <- "/Users/ricard/10x_gastrulation_DNMTs/dimensionality_reduction/dimensionality_reduction_sce.R"
} else {
  source("/homes/ricard/10x_gastrulation_DNMTs/settings.R")
  io$script <- "/homes/ricard/10x_gastrulation_DNMTs/dimensionality_reduction/dimensionality_reduction_sce.R"
  io$tmpdir <- paste0(io$basedir,"/results/dimensionality_reduction/tmp"); dir.create(io$tmpdir)
}
io$outdir <- paste0(io$basedir,"/results/dimensionality_reduction")

####################
## Define options ##
####################

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
  # "15_E8_5_D3A_WT_D3B_WT_L007",
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
  # "SIGAH9_Dnmt3a_KO_Dnmt3b_Het_L001",
  # "SIGAG5_9_dnmt3ab_DKO_L005"
)

opts$batches <- "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001"

# Test mode (subsetting cells)?
opts$test_mode <- FALSE


# Number of highly variable genes
opts$features <- 1000

# Number of PCs
opts$npcs <- 30

#########
## Run ##
#########

for (i in opts$batches) {
  
  # Define LSF command
  if (grepl("ricard",Sys.info()['nodename'])) {
    lsf <- ""
  } else if (grepl("ebi",Sys.info()['nodename'])) {
    lsf <- sprintf("bsub -M 10000 -n 1 -o %s/%s.txt", io$tmpdir,paste(i,collapse=" "))
  }
  cmd <- sprintf("%s Rscript %s --batches %s --features %d --npcs %d --outdir %s", lsf, io$script, i, opts$features, opts$npcs, io$outdir)
  
  if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test")
  
  # Run
  print(cmd)
  system(cmd)
}
