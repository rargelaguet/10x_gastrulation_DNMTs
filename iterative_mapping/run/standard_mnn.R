# Load libraries
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(argparse))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--atlas_stages',    type="character",   nargs='+',  help='Atlas stage(s)')
p$add_argument('--query_batches',   type="character",   nargs='+',  help='Query batch(es)')
p$add_argument('--test',            action = "store_true",  help = 'Testing mode')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$atlas_stages <- c(
#   # "E6.5"
#   # "E6.75",
#   # "E7.0",
#   # "E7.25",
#   # "E7.5",
#   # "E7.75",
#   # "E8.0",
#   # "E8.25",
#   "E8.5"
#   # "mixed_gastrulation"
# )
# args$query_batches <- c(
#   "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001"
#   # "SIGAB6_E85_3_Dnmt3aWT_Dnmt3b_WT_L002",
#   # "SIGAC6_E85_5_Dnmt3aKO_Dnmt3b_Het_L003",
#   # "SIGAD6_E85_8_Dnmt3aHet_Dnmt3b_KO_L004",
#   # "15_E8_5_D3A_WT_D3B_WT_L007",
#   # "17_E8_5_D3A_KO_D3B_WT_L008",
#   # "2_E8_5_D3A_WT_D3B_KO_L003",
#   # "3_E8_5_D3A_HET_D3B_WT_L004",
#   # "7_E8_5_D3A_WT_D3B_KO_L005",
#   # "8_E8_5_D3A_KO_D3B_KO_L006"
# )
# args$test <- TRUE
## END TEST ##

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
  source("/Users/ricard/10x_gastrulation_DNMTs/iterative_mapping/run/utils.R")
  io$atlas.marker_genes <- "/Users/ricard/data/gastrulation10x/results/marker_genes/E8.5/marker_genes.txt.gz"
  io$script_load_data <- "/Users/ricard/10x_gastrulation_DNMTs/iterative_mapping/run/load_data.R"
} else {
  source("/homes/ricard/10x_gastrulation_DNMTs/settings.R")
  source("/homes/ricard/10x_gastrulation_DNMTs/iterative_mapping/run/utils.R")
  io$atlas.marker_genes <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/results/marker_genes/E8.5/marker_genes.txt.gz"
  io$script_load_data <- "/homes/ricard/10x_gastrulation_DNMTs/iterative_mapping/run/load_data.R"
}
io$path2atlas <- io$atlas.basedir
io$path2query <- io$basedir
io$outdir <- paste0(io$basedir,"/results/iterative_mapping")

if (isTRUE(args$test)) print("Test mode activated...")

###############
## Load data ##
###############

source(io$script_load_data)

#############
## Run MNN ##
#############

# Joint normalisation
sce.all <- joint.normalisation(sce.query, sce.atlas, cosineNorm = FALSE)

# Select HVGs
# genes <- getHVGs(sce.atlas, block=as.factor(sce.atlas$sample), p.value = 0.10)
genes <- getHVGs(sce.all, block=sce.all$block, p.value = 0.05)

# MNN
mapping.dt <- mnn.fn(sce.all, sce.query, sce.atlas, genes = genes, npcs = 50, k = 25)

# table(mapping.dt$celltype_mapped)
# foo <- mapping.dt %>% merge(meta_query[,c("cell","celltype.mapped","celltype.score")] %>% setnames(c("cell","celltype_old","score_old")))

##########
## Save ##
##########

fwrite(mapping.dt, sprintf("%s/%s_standard_mnn.txt.gz",io$outdir,paste(args$query_batches,collapse="_")), sep="\t")

