suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scds))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(argparse))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--samples',                 type="character",   nargs='+',     help='sample(es)')
p$add_argument('--hybrid_score_threshold',  type="double",      default=1.0,   help='Hybrid score threshold')
p$add_argument('--test',                    action = "store_true",             help='Testing mode')
p$add_argument('--outdir',                  type="character",                  help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args <- list()
# args$samples <- c(
#   "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001"
# )
# 
# args$hybrid_score_threshold <- 1.0
# args$npcs <- 30
# args$test <- TRUE
## END TEST ##

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
} else {
  source("/homes/ricard/10x_gastrulation_DNMTs/settings.R")
}

if (isTRUE(args$test)) print("Test mode activated...")

############################
## Update sample metadata ##
############################

sample_metadata <- fread(io$metadata) %>%
  .[pass_QC==TRUE & sample%in%args$samples]
table(sample_metadata$sample)

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- load_SingleCellExperiment(args$sce, cells=sample_metadata$cell, normalise = TRUE)
dim(sce)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

#############################
## Calculate doublet score ##
#############################

sce <- cxds_bcds_hybrid(sce, estNdbl=TRUE)

dt <- colData(sce) %>%
    .[,c("sample","cxds_score", "cxds_call", "bcds_score", "bcds_call", "hybrid_score", "hybrid_call")] %>%
    as.data.frame %>% tibble::rownames_to_column("cell") %>% as.data.table

# Call doublets
dt[,hybrid_call:=hybrid_score>args$hybrid_score_threshold]
table(dt$hybrid_call)

##########
## Save ##
##########

io$outfile <- sprintf("%s/%s_%s.txt.gz",args$outdir, paste(args$samples,collapse="-"),round(args$hybrid_score_threshold,2))
fwrite(dt, io$outfile, sep="\t", na="NA", quote=F)
