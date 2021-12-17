#####################
## Define settings ##
#####################

io <- list(); opts <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$standard.mnn.script <- "/Users/ricard/10x_gastrulation_DNMTs/iterative_mapping/run/standard_mnn.R"
  io$iterative.mnn.script <- "/Users/ricard/10x_gastrulation_DNMTs/iterative_mapping/run/iterative_mnn.R"
} else if(grepl("ebi",Sys.info()['nodename'])){
  io$standard.mnn.script <- "/homes/ricard/10x_gastrulation_DNMTs/iterative_mapping/run/standard_mnn.R"
  io$iterative.mnn.script <- "/homes/ricard/10x_gastrulation_DNMTs/iterative_mapping/run/iterative_mnn.R"
  io$tmpdir <- "/hps/nobackup2/research/stegle/users/ricard/10x_gastrulation_DNMTs/results/iterative_mapping/tmp"
} 

# Atlas stages
opts$atlas_stages <- c(
  # "E6.5",
  # "E6.75",
  # "E7.0",
  # "E7.25",
  "E7.5",
  "E7.75",
  "E8.0",
  "E8.25",
  "E8.5"
  # "mixed_gastrulation"
)

# Query batches
opts$query_batches <- c(
  # "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001",
  # "SIGAB6_E85_3_Dnmt3aWT_Dnmt3b_WT_L002",
  # "SIGAC6_E85_5_Dnmt3aKO_Dnmt3b_Het_L003",
  # "SIGAD6_E85_8_Dnmt3aHet_Dnmt3b_KO_L004",
  "15_E8_5_D3A_WT_D3B_WT_L007",
  # "17_E8_5_D3A_KO_D3B_WT_L008",
  # "2_E8_5_D3A_WT_D3B_KO_L003",
  # "3_E8_5_D3A_HET_D3B_WT_L004",
  # "7_E8_5_D3A_WT_D3B_KO_L005",
  "8_E8_5_D3A_KO_D3B_KO_L006"
)

# Test mode (subset cells)?
opts$test <- FALSE


#########
## Run ##
#########

for (i in opts$query_batches) {
  # LSF
  if (grepl("ricard",Sys.info()['nodename'])) {
    lsf <- ""
  } else {
    lsf <- sprintf("bsub -M 60000 -n 1 -q research-rh74 -o %s/%s.txt", io$tmpdir, i)
  }

  # Run standard MNN
  cmd <- sprintf("%s Rscript %s --query_batches %s --atlas_stages %s", lsf, io$standard.mnn.script, i, paste(opts$atlas_stages, collapse=" "))
  if (isTRUE(opts$test)) cmd <- paste0(cmd, " --test")
  system(cmd)

  # Run tree-guided MNN
  cmd <- sprintf("%s Rscript %s --query_batches %s --atlas_stages %s", lsf, io$iterative.mnn.script, i, paste(opts$atlas_stages, collapse=" "))
  if (isTRUE(opts$test)) cmd <- paste0(cmd, " --test")
  system(cmd)
}
