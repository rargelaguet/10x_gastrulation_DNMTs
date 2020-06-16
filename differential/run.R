#########
## I/O ##
#########

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
  io$script <- "/Users/ricard/10x_gastrulation_DNMTs/differential/differential.R"
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/10x_gastrulation_DNMTs/settings.R")
  io$script <- "/homes/ricard/10x_gastrulation_DNMTs/differential/differential.R"
  io$tmpdir <- "/hps/nobackup2/research/stegle/users/ricard/10x_gastrulation_DNMTs/results/differential/tmp"; dir.create(io$tmpdir, showWarnings=F)
} else {
  stop("Computer not recognised")
}
io$outdir <- paste0(io$basedir,"/results/differential"); dir.create(io$outdir, showWarnings = F)

#############
## Options ##
#############

# Testing mode
opts$test_mode <- FALSE

# Define groups
# opts$query.groups <- c(
#   "E8.5_Dnmt3aKO_Dnmt3bWT", 
#   "E8.5_Dnmt3aWT_Dnmt3bWT", 
#   "E8.5_Dnmt3aKO_Dnmt3bHET", 
#   "E8.5_Dnmt3aHET_Dnmt3bKO",
#   "E12.5_Dnmt3aWT_Dnmt3bHET", 
#   "E12.5_Dnmt3aWT_Dnmt3bKO" 
# )

opts$query.groups <- opts$classes[opts$classes!="E8.5_Dnmt3aWT_Dnmt3bWT"] %>% head(n=1)
opts$reference.groups <- c("E8.5_Dnmt3aWT_Dnmt3bWT")
# opts$groupA <- "E8.5_Dnmt3aWT_Dnmt3bWT"
# opts$groupB <- "E85_Dnmt3aWT_Dnmt3bWT"


  # Define cell types with sufficient number of cells for differential comparison
opts$min.cells <- 25

#########
## Run ##
#########

i <- opts$query.groups[1]
for (i in opts$query.groups) {
  
  # Define groups
  opts$groupA <- i
  opts$groupB <- opts$reference.groups
  
  # Extract cell types
  opts$celltypes <- sample_metadata %>% 
    .[!is.na(celltype.mapped2) & class%in%c(opts$groupA,opts$groupB)] %>%
    .[,.N,by=c("class","celltype.mapped2")] %>%
    .[N>=opts$min.cells] %>% .[,.N,by="celltype.mapped2"] %>% .[N==2,celltype.mapped2]
  
    j <- opts$celltypes[1]
    for (j in opts$celltypes) {
      outfile <- sprintf("%s/%s_vs_%s_%s.txt.gz", io$outdir,opts$groupA,opts$groupB,j)

      # Define LSF command
      if (grepl("ricard",Sys.info()['nodename'])) {
        lsf <- ""
      } else if (grepl("ebi",Sys.info()['nodename'])) {
        lsf <- sprintf("bsub -M 15000 -n 1 -q research-rh74 -o %s/%s_vs_%s_%s.txt", io$tmpdir,opts$groupA,opts$groupB,j)
      }
      cmd <- sprintf("%s Rscript %s --groupA %s --groupB %s --celltype %s --outfile %s", lsf, io$script, opts$groupA, opts$groupB, j, outfile)
      if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test_mode")

      # Run
      print(cmd)
      system(cmd)
    }
}

