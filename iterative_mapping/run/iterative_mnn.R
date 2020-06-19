suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(edgeR))

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
# args$test <- FALSE
## END TEST ##

#####################
## Define settings ##
#####################

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

########################################################
## Define distance matrix for hierarchical clustering ##
########################################################

opts$celltypes <- which(table(meta_atlas$celltype)>25) %>% names %>% stringr::str_replace_all("_", " ")
# opts$celltypes <- unique(sample_metadata_atlas$celltype) 

dist <- fread(paste0(io$atlas.basedir,"/results/phylogenetic_tree/PAGA_distances.csv.gz")) %>%
  as.data.frame %>% tibble::column_to_rownames("V1") %>% as.matrix %>%
  .[opts$celltypes,opts$celltypes] %>%
  as.dist

#######################
## Recursive mapping ##
#######################

sce.query$celltype_mapped <- paste(opts$celltypes,collapse="%")

while (any(grepl("%",sce.query$celltype_mapped))) {
  print(table(sce.query$celltype_mapped))
  mapping.dt <- recursive.fn(sce.query, sce.atlas, dist, cosineNorm = TRUE)
  ids <- match(mapping.dt$cell,colnames(sce.query))
  sce.query$celltype_mapped[ids] <- mapping.dt$celltype_mapped
  sce.query$celltype_score[ids] <- mapping.dt$celltype_score
  sce.query$stage_mapped[ids] <- mapping.dt$stage_mapped
  sce.query$stage_score[ids] <- mapping.dt$stage_score
}

##########
## Save ##
##########

mapping.dt <- data.table(
  cell = colnames(sce.query), 
  celltype_mapped = sce.query$celltype_mapped,
  celltype_score = sce.query$celltype_score,
  stage_mapped = sce.query$stage_mapped,
  stage_score = sce.query$stage_score
)
# foo <- mapping.dt %>% merge(meta_query[,c("cell","celltype.mapped","celltype.score")] %>% setnames(c("cell","celltype_old","score_old")))

fwrite(mapping.dt, sprintf("%s/%s_iterative_mnn.txt.gz",io$outdir,paste(args$test_samples,collapse="_")), sep="\t")
