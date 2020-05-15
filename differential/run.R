#########
## I/O ##
#########

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
  io$script <- "/Users/ricard/10x_gastrulation_DNMTs/differential/differential.R"
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/10x_gastrulation_DNMTs/settings.R")
  io$script <- "/homes/ricard/10x_gastrulation_DNMTs/differential/differential.R"
  io$tmpdir <- "/hps/nobackup2/research/stegle/users/ricard/10x_gastrulation_DNMTs/results/second_batch/differential/tmp"; dir.create(io$tmpdir, showWarnings=F)
} else {
  stop("Computer not recognised")
}
io$outdir <- paste0(io$basedir,"/results/differential"); dir.create(io$outdir, showWarnings = F)

#############
## Options ##
#############

# Statistical test
opts$statistical.test <- "edgeR"

# Testing mode
opts$test_mode <- TRUE

# Define groups
opts$groupA <- "Dnmt3aWT_Dnmt3bHET"
opts$groupB <- "Dnmt3aWT_Dnmt3bKO"

# Define cell types with sufficient number of cells for differential comparison
opts$min.cells <- 25
opts$celltypes <- sample_metadata %>%
  .[!is.na(celltype.mapped) & class%in%c(opts$groupA,opts$groupB),.N,by=c("class","celltype.mapped")] %>%
  .[N>=opts$min.cells] %>% .[,.N,by="celltype.mapped"] %>% .[N==2,celltype.mapped]

#########
## Run ##
#########

for (i in head(opts$celltypes,n=3)) {
# for (i in opts$celltypes.1) {
    outfile <- sprintf("%s/%s_WT_vs_TET_TKO.txt.gz", io$outdir,i)
    
    # Define LSF command
    if (grepl("ricard",Sys.info()['nodename'])) {
      lsf <- ""
    } else if (grepl("ebi",Sys.info()['nodename'])) {
      lsf <- sprintf("bsub -M 15000 -n 1 -q research-rh74 -o %s/%s_vs_%s.txt", io$tmpdir,opts$groupA,opts$groupB)
    }
    cmd <- sprintf("%s Rscript %s --groupA %s --groupB %s --celltype %s --test %s --outfile %s", lsf, io$script, opts$groupA, opts$groupB, i, opts$statistical.test, outfile)
    if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test_mode")
    
    # Run
    print(cmd)
    system(cmd)
}

