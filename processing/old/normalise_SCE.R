suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scran))

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/10x_gastrulation_DNMTs/settings.R")
}

## I/O ##

io$outdir <- paste0(io$basedir,"/results/individual_genes")
