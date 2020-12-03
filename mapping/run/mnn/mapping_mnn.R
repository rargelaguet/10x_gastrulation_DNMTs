suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(argparse))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--atlas_stages',    type="character",   nargs='+',  help='Atlas stage(s)')
p$add_argument('--query_batches',   type="character",   nargs='+',  help='Query batch(es)')
p$add_argument('--test',            action = "store_true",          help='Testing mode')
p$add_argument('--outdir',          type="character",               help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$atlas_stages <- c(
#   "E6.5",
#   "E6.75",
#   "E7.0",
#   "E7.25",
#   "E7.5",
#   "E7.75",
#   "E8.0",
#   "E8.25",
#   "E8.5",
#   "mixed_gastrulation"
# )
# 
# args$query_batches <- c(
#   "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001"
# )
# 
# args$test <- TRUE
## END TEST ##

################
## Define I/O ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
  source("/Users/ricard/10x_gastrulation_DNMTs/mapping/run/mapping_functions.R")
} else {
  source("/homes/ricard/10x_gastrulation_DNMTs/settings.R")
  source("/homes/ricard/10x_gastrulation_DNMTs/mapping/run/mapping_functions.R")
}
io$path2atlas <- io$atlas.basedir
io$path2query <- io$basedir

if (isTRUE(args$test)) print("Test mode activated...")

####################
## Define options ##
####################

opts$npcs <- 50
opts$k <- 25

################
## Load atlas ##
################

# Load SingleCellExperiment
sce_atlas  <- readRDS(io$atlas.sce)

# Load cell metadata
meta_atlas <- fread(io$atlas.metadata) %>%
  .[stripped==F & doublet==F & stage%in%args$atlas_stages]

# Filter
if (isTRUE(args$test)) meta_atlas <- head(meta_atlas,n=1000)
sce_atlas <- sce_atlas[,meta_atlas$cell] 

################
## Load query ##
################

# Load SingleCellExperiment
io$sce <- "/hps/nobackup2/research/stegle/users/ricard/10x_gastrulation_DNMTs/processed/sixth_batch/SingleCellExperiment.rds"
sce_query <- readRDS(io$sce)

# Load cell metadata
io$metadata <- "/hps/nobackup2/research/stegle/users/ricard/10x_gastrulation_DNMTs/processed/sixth_batch/metadata.txt.gz"
meta_query <- fread(io$metadata) %>% .[pass_QC==T & batch%in%args$query_batches]

# Filter
if (isTRUE(args$test)) meta_query <- head(meta_query,n=1000)
sce_query <- sce_query[,meta_query$cell]

#############
## Prepare ## 
#############

# Filter out non-expressed genes
sce_query <- sce_query[rowSums(counts(sce_query))>10,]
sce_atlas <- sce_atlas[rowSums(counts(sce_atlas))>10,]

# Rename ensemble IDs to gene names
gene_metadata <- fread(io$gene_metadata) %>% .[,c("ens_id","symbol")] %>%
  .[symbol!="" & ens_id%in%rownames(sce_atlas)] %>%
  .[!duplicated(symbol)]

sce_atlas <- sce_atlas[rownames(sce_atlas)%in%gene_metadata$ens_id,]
foo <- gene_metadata$symbol; names(foo) <- gene_metadata$ens_id
rownames(sce_atlas) <- foo[rownames(sce_atlas)]

# Sanity cehcks
stopifnot(sum(is.na(rownames(sce_atlas)))==0)
stopifnot(sum(duplicated(rownames(sce_atlas)))==0)

# Intersect genes
genes.intersect <- intersect(rownames(sce_query), rownames(sce_atlas))

# Remove mitochondrial genes
genes.intersect <- genes.intersect[grep("mt-",genes.intersect,invert = T)]
genes.intersect <- genes.intersect[grep("Rik",genes.intersect,invert = T)]

# Subset SingleCellExperiment objects
sce_query  <- sce_query[genes.intersect,]
sce_atlas <- sce_atlas[genes.intersect,]

#########
## Map ##
#########

mapping  <- mapWrap(
  atlas_sce = sce_atlas, atlas_meta = meta_atlas,
  map_sce = sce_query, map_meta = meta_query, 
  # genes = marker_genes, 
  npcs = opts$npcs, k = opts$k
)

##########
## Save ##
##########

outfile <- sprintf("%s/mapping_mnn_%s.rds",args$outdir,paste(args$query_batches,collapse="-"))
saveRDS(mapping, outfile)
