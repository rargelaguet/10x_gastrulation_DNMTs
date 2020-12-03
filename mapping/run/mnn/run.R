
################
## Define I/O ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
  io$script <- "/Users/ricard/10x_gastrulation_DNMTs/mapping/run/mnn/mapping_mnn.R"
} else {
  source("/homes/ricard/10x_gastrulation_DNMTs/settings.R")
  io$script <- "/homes/ricard/10x_gastrulation_DNMTs/mapping/run/mnn/mapping_mnn.R"
  io$tmpdir <- paste0(io$basedir,"/results/mapping/tmp"); dir.create(io$tmpdir)
}
io$outdir <- paste0(io$basedir,"/results/mapping")

####################
## Define options ##
####################

opts$atlas_stages <- c(
  "E6.5",
  "E6.75",
  "E7.0",
  "E7.25",
  "E7.5",
  "E7.75",
  "E8.0",
  "E8.25",
  "E8.5",
  "mixed_gastrulation"
)

opts$batches <- c(
  # E12.5  
  # "E125_DNMT3A_HET_A_L001",
  # "E125_DNMT3A_HET_A_L003",
  # "E125_DNMT3A_KO_B_L002",
  # "E125_DNMT3A_KO_E_L004",
  # "A_E12_5_D3a_Het_L001",
  # "B_E12_5_D3a_KO_L002",
  
  
  # E8.5  
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
  "SIGAH9_Dnmt3a_KO_Dnmt3b_Het_L001"
)

# opts$batches <- "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001"

# Test mode (subsetting cells)?
opts$test_mode <- FALSE

if (opts$test_mode) {
  opts$memory <- 25000
} else {
  opts$memory <- 40000
}
#########
## Run ##
#########

for (i in opts$batches) {
  
  # Define LSF command
  if (grepl("ricard",Sys.info()['nodename'])) {
    lsf <- ""
  } else if (grepl("ebi",Sys.info()['nodename'])) {
    lsf <- sprintf("bsub -M %d -n 1 -o %s/%s.txt", opts$memory, io$tmpdir,paste(i,collapse=" "))
  }
  cmd <- sprintf("%s Rscript %s --atlas_stages %s --query_batches %s --outdir %s", lsf, io$script, paste(opts$atlas_stages,collapse=" "), paste(i,collapse=" "),io$outdir)
  if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test")
  
  # Run
  print(cmd)
  system(cmd)
}
