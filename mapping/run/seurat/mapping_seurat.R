suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(argparse))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--query_sce',       type="character",               help='SingleCellExperiment file for the query')
p$add_argument('--atlas_sce',       type="character",               help='SingleCellExperiment file for the atlas')
p$add_argument('--query_metadata',  type="character",               help='metadata file for the query')
p$add_argument('--atlas_metadata',  type="character",               help='metadata file for the atlas')
p$add_argument('--atlas_stages',    type="character",   nargs='+',  help='Atlas stage(s)')
p$add_argument('--query_samples',   type="character",   nargs='+',  help='Query batch(es)')
p$add_argument('--npcs',            type="integer",                 help='Number of principal components')
p$add_argument('--n_neighbours',    type="integer",                 help='Number of neighbours')
p$add_argument('--use_marker_genes',action = "store_true",          help='Use marker genes?')
p$add_argument('--test',            action = "store_true",          help='Testing mode')
p$add_argument('--outfile',          type="character",               help='Output file')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

# I/O
io$path2atlas <- io$atlas.basedir
io$path2query <- io$basedir

# Load mapping functions
# source(here::here("rna/mapping/run/mapping_functions.R"))


## START TEST ##
# args$query_sce <- io$rna.sce
# args$atlas_sce <- io$rna.atlas.sce
# args$query_metadata <- io$metadata
# args$atlas_metadata <- io$rna.atlas.metadata
# args$query_samples <- opts$samples[1:2]
# args$atlas_stages <- c("E6.5")
# args$npcs <- 50
# args$n_neighbours <- 25
# args$use_marker_genes <- FALSE
# args$test <- TRUE
# args$outfile <- paste0(io$basedir,"/results/rna/mapping/test")
## END TEST ##

if (isTRUE(args$test)) print("Test mode activated...")

################
## Load query ##
################

# Load cell metadata
meta_query <- fread(args$query_metadata) %>% .[pass_rnaQC==T & sample%in%args$query_samples]
if (isTRUE(args$test)) meta_query <- head(meta_query,n=1000)

# Load SingleCellExperiment
sce_query <- load_SingleCellExperiment(args$query_sce, cells = meta_query$cell, remove_non_expressed_genes = TRUE)

# Update colData
tmp <- meta_query %>% .[cell%in%colnames(sce_query)] %>% setkey(cell) %>% .[colnames(sce_query)]
stopifnot(tmp$cell == colnames(sce_query))
colData(sce_query) <- tmp %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
  .[colnames(sce_query),] %>% DataFrame()

################
## Load atlas ##
################

# Load cell metadata
meta_atlas <- fread(args$atlas_metadata) %>%
  .[stripped==F & doublet==F & stage%in%args$atlas_stages]

# Filter
if (isTRUE(args$test)) meta_atlas <- head(meta_atlas,n=1000)

# Load SingleCellExperiment
sce_atlas <- load_SingleCellExperiment(args$atlas_sce, normalise = TRUE, cells = meta_atlas$cell, remove_non_expressed_genes = TRUE)

# Update colData
tmp <- meta_atlas %>% .[cell%in%colnames(sce_atlas)] %>% setkey(cell) %>% .[colnames(sce_atlas)]
stopifnot(tmp$cell == colnames(sce_atlas))
colData(sce_atlas) <- tmp %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
  .[colnames(sce_atlas),] %>% DataFrame()

#############
## Prepare ## 
#############

# Rename ensemble IDs to gene names in the atlas
gene_metadata <- fread(io$gene_metadata) %>% .[,c("chr","ens_id","symbol")] %>%
  .[symbol!="" & ens_id%in%rownames(sce_atlas)] %>%
  .[!duplicated(symbol)]

sce_atlas <- sce_atlas[rownames(sce_atlas)%in%gene_metadata$ens_id,]
foo <- gene_metadata$symbol; names(foo) <- gene_metadata$ens_id
rownames(sce_atlas) <- foo[rownames(sce_atlas)]

# Sanity checks
stopifnot(sum(is.na(rownames(sce_atlas)))==0)
stopifnot(sum(duplicated(rownames(sce_atlas)))==0)

#####################
## Define gene set ##
#####################

# Intersect genes
genes.intersect <- intersect(rownames(sce_query), rownames(sce_atlas))

# Filter some genes manually
genes.intersect <- genes.intersect[grep("^Rik|Rik$|^mt-|^Rps-|^Rpl-|^Gm",genes.intersect,invert=T)]
genes.intersect <- genes.intersect[!genes.intersect=="Xist"]
genes.intersect <- genes.intersect[!genes.intersect%in%gene_metadata[chr=="chrY",symbol]]

# Subset marker genes
if (args$use_marker_genes) {
  marker_genes.dt <- fread(io$rna.atlas.marker_genes)
  genes.intersect <- genes.intersect[genes.intersect%in%unique(marker_genes.dt$gene)]
}

# Subset Seurat objects
sce_query  <- sce_query[genes.intersect,]
sce_atlas <- sce_atlas[genes.intersect,]

###########################
## Create Seurat objects ##
###########################

seurat_atlas   <- CreateSeuratObject(as.matrix(counts(sce_atlas)), project = "ATLAS", min.cells = 25)
seurat_atlas@meta.data$map <- "ATLAS"
seurat_atlas   <- NormalizeData(seurat_atlas)
rownames(meta_atlas) <- meta_atlas$cell
seurat_atlas   <- AddMetaData(object = seurat_atlas, metadata = meta_atlas)
# seurat_atlas   <- ScaleData(seurat_atlas, display.progress = TRUE, vars.to.regress = "sample")
seurat_atlas   <- ScaleData(seurat_atlas)

seurat_query <- CreateSeuratObject(as.matrix(counts(sce_query)), project = "QUERY", min.cells = 25)
seurat_query@meta.data$map <- "QUERY"
seurat_query <- NormalizeData(seurat_query)
seurat_query <- ScaleData(seurat_query)

#########
## Map ##
#########

# Find anchors
if (args$use_marker_genes) {
  anchors <- FindTransferAnchors(reference = seurat_atlas, query = seurat_query, features=rownames(seurat_query), dims = 1:args$npcs, k.anchor = args$n_neighbours)
} else {
  seurat_atlas <- FindVariableFeatures(seurat_atlas, selection.method = "vst")
  seurat_query <- FindVariableFeatures(seurat_query, selection.method = "vst")
  anchors <- FindTransferAnchors(reference = seurat_atlas, query = seurat_query, dims = 1:args$npcs, k.anchor = args$n_neighbours)
}

# Transfer data
predictions <- TransferData(anchorset = anchors, refdata = seurat_atlas$celltype, dims = 1:args$npcs, k.weight = args$n_neighbours)
seurat_query <- AddMetaData(object = seurat_query, metadata = predictions)

##########
## Save ##
##########

mapping.dt <- predictions %>% as.data.table(keep.rownames = T) %>%
  .[,c("rn","predicted.id","prediction.score.max")] %>%
    setnames(c("cell","celltype.mapped","celltype.score")) %>%
  .[,celltype.score:=round(celltype.score,2)]

# outfile <- sprintf("%s/mapping_seurat_%s.txt.gz",args$outdir,paste(args$query_samples,collapse="-"))
# fwrite(mapping.dt, outfile, sep="\t")
fwrite(mapping.dt, args$outfile, sep="\t")
