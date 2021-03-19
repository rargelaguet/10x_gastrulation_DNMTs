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
p$add_argument('--query_sce',       type="character",               help='SingleCellExperiment file for the query')
p$add_argument('--atlas_sce',       type="character",               help='SingleCellExperiment file for the atlas')
p$add_argument('--query_metadata',  type="character",               help='metadata file for the query')
p$add_argument('--atlas_metadata',  type="character",               help='metadata file for the atlas')
p$add_argument('--npcs',            type="integer",                 help='Number of principal components')
p$add_argument('--n_neighbours',    type="integer",                 help='Number of neighbours')
p$add_argument('--test',            action = "store_true",          help='Testing mode')
p$add_argument('--outdir',          type="character",               help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
  source("/Users/ricard/10x_gastrulation_DNMTs/mapping/run/mapping_functions.R")
} else {
  source("/homes/ricard/10x_gastrulation_DNMTs/settings.R")
  source("/homes/ricard/10x_gastrulation_DNMTs/mapping/run/mapping_functions.R")
}
io$path2atlas <- io$atlas.basedir
io$path2query <- io$basedir

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
# args$query_sce <- io$sce
# args$atlas_sce <- io$atlas.sce
# args$query_metadata <- io$metadata
# args$atlas_metadata <- io$atlas.metadata
# args$test <- TRUE
# args$npcs <- 50
# args$n_neighbours <- 25
## END TEST ##

if (isTRUE(args$test)) print("Test mode activated...")

####################
## Define options ##
####################

################
## Load atlas ##
################

# Load cell metadata
meta_atlas <- fread(args$atlas_metadata) %>%
  .[stripped==FALSE & doublet==FALSE & stage%in%args$atlas_stages]

# Filter
if (isTRUE(args$test)) meta_atlas <- head(meta_atlas,n=1000)

# Load SingleCellExperiment
sce_atlas <- load_SingleCellExperiment(args$atlas_sce, cells=meta_atlas$cell, normalise = TRUE, remove_non_expressed_genes = TRUE)

################
## Load query ##
################


# Load cell metadata
meta_query <- fread(args$query_metadata) %>% .[pass_QC==TRUE & batch%in%args$query_batches]

# Filter
if (isTRUE(args$test)) meta_query <- head(meta_query,n=1000)

# Load SingleCellExperiment
sce_query <- load_SingleCellExperiment(args$query_sce, cells = meta_query$cell, normalise = FALSE, remove_non_expressed_genes = TRUE)

#############
## Prepare ## 
#############

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
# genes.intersect <- genes.intersect[grep("Rik",genes.intersect,invert = T)]

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
  npcs = args$npcs, 
  k = args$n_neighbours
)

##########
## Save ##
##########

outfile <- sprintf("%s/mapping_mnn_%s.rds",args$outdir,paste(args$query_batches,collapse="-"))
saveRDS(mapping, outfile)

##########
## TEST ##
##########

# foo <- mapping$mapping %>% as.data.table() %>% .[,c("cell","celltype.mapped","celltype.score")]
# bar <- fread("/Users/ricard/data/gastrulation_multiome_10x/results/rna/mapping/sample_metadata_after_mapping.txt.gz") %>% 
#   .[,c("cell","celltype.mapped","celltype.score")]
# foobar <- merge(foo,bar,by=c("cell"))
# foobar[celltype.mapped.x!=celltype.mapped.y] %>% View
# mean(foobar$celltype.mapped.x==foobar$celltype.mapped.y)
# fwrite(foo,"/Users/ricard/data/gastrulation_multiome_10x/results/rna/mapping/E7.5_rep2_test.txt.gz")
