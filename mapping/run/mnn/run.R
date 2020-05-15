
################
## Define I/O ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
  io$script <- "/Users/ricard/10x_gastrulation_DNMTs/mapping/run/mnn/mapping_mnn.R"
} else {
  source("/homes/ricard/10x_gastrulation_DNMTs/settings.R")
  io$script <- "/homes/ricard/10x_gastrulation_DNMTs/mapping/run/mnn/mapping_mnn.R"
  io$tmpdir <- paste0(io$basedir,"/results/second_batch/mapping/tmp"); dir.create(io$tmpdir)
}
io$outdir <- paste0(io$basedir,"/results/second_batch/mapping")

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

opts$query_batches <- c(
  "E75_TET_TKO_L002", 
  "E75_WT_Host_L001", 
  "E85_Rep1_TET_TKO_L004", 
  "E85_Rep2_TET_TKO_L006", 
  "E85_Rep1_WT_Host_L003",
  "E85_Rep2_WT_Host_L005"
  # "E125_DNMT3A_HET_A_L001", 
  # "E125_DNMT3A_HET_A_L003",
  # "E125_DNMT3A_KO_B_L002", 
  # "E125_DNMT3A_KO_E_L004"
)

# Test mode (subsetting cells)?
opts$test_mode <- FALSE

#########
## Run ##
#########

for (i in opts$query_batches) {
  
  # Define LSF command
  if (grepl("ricard",Sys.info()['nodename'])) {
    lsf <- ""
  } else if (grepl("ebi",Sys.info()['nodename'])) {
    lsf <- sprintf("bsub -M 70000 -n 1 -o %s/%s.txt", io$tmpdir,paste(i,collapse=" "))
  }
  cmd <- sprintf("%s Rscript %s --atlas_stages %s --query_batches %s --outdir %s", lsf, io$script, paste(opts$atlas_stages,collapse=" "), paste(i,collapse=" "),io$outdir)
  if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test")
  
  # Run
  print(cmd)
  system(cmd)
}
