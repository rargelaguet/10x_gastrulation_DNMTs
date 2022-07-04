here::i_am("mapping/run/mnn/mapping_mnn.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(batchelor))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--atlas_stages',    type="character",   nargs='+',  help='Atlas stage(s)')
p$add_argument('--query_samples',   type="character",   nargs='+',  help='Query batch(es)')
p$add_argument('--query_sce',       type="character",               help='SingleCellExperiment file for the query')
p$add_argument('--atlas_sce',       type="character",               help='SingleCellExperiment file for the atlas')
p$add_argument('--query_metadata',  type="character",               help='metadata file for the query')
p$add_argument('--atlas_metadata',  type="character",               help='metadata file for the atlas')
p$add_argument('--npcs',            type="integer",                 help='Number of principal components')
p$add_argument('--n_neighbours',    type="integer",                 help='Number of neighbours')
p$add_argument('--use_marker_genes',action = "store_true",          help='Use marker genes?')
p$add_argument('--cosine_normalisation',      action = "store_true",          help='Use cosine normalisation?')
p$add_argument('--test',            action = "store_true",          help='Testing mode')
p$add_argument('--outfile',          type="character",               help='Output file')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

# I/O
io$path2atlas <- io$atlas.basedir
io$path2query <- io$basedir

# Load mapping functions
source(here::here("mapping/run/mnn/mapping_functions.R"))

## START TEST ##
args$atlas_stages <- c("E6.5","E6.75","E7.0","E7.25","E7.5","E7.75","E8.0","E8.25","E8.5")
args$query_samples <- "SIGAH10_Dnmt3ab_WT_L002"# opts$samples[1]
args$query_sce <- paste0(io$basedir,"/processed/SingleCellExperiment.rds")
args$atlas_sce <- io$atlas.sce
args$query_metadata <- paste0(io$basedir,"/results/qc/sample_metadata_after_qc.txt.gz")
args$atlas_metadata <- io$atlas.metadata
args$test <- FALSE
args$npcs <- 50
args$n_neighbours <- 25
args$use_marker_genes <- FALSE
args$cosine_normalisation <- FALSE
args$outdir <- paste0(io$basedir,"/results/mapping/test")
## END TEST ##

if (isTRUE(args$test)) print("Test mode activated...")

################
## Load query ##
################

# Load cell metadata
meta_query <- fread(args$query_metadata) %>% 
  # .[pass_rnaQC==TRUE & doublet_call==FALSE & sample%in%args$query_samples]
  .[pass_rnaQC==TRUE & sample%in%args$query_samples]
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
  .[stripped==F & doublet==F & stage%in%args$atlas_stages] %>%
  .[,sample:=factor(sample)]

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

# Sanity cehcks
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

# Subset SingleCellExperiment objects
sce_query  <- sce_query[genes.intersect,]
sce_atlas <- sce_atlas[genes.intersect,]

#######################
## Feature selection ##
#######################

if (args$use_marker_genes) {
  # Load marker genes
  marker_genes.dt <- fread(io$rna.atlas.marker_genes)
  genes_to_use <- genes.intersect[genes.intersect%in%unique(marker_genes.dt$gene)]
} else {
  # Load gene statistics from the atlas
  gene_stats.dt <- fread(paste0(io$atlas.basedir,"/results/gene_statistics/gene_statistics.txt.gz")) %>%
    .[gene%in%genes.intersect]
  genes_to_use <- gene_stats.dt %>% setorder(-var_pseudobulk, na.last = T) %>% head(n=2500) %>% .$gene  
  
  # Calculate mean-variance relationship and extract HVGs
  # decomp <- modelGeneVar(sce_atlas, block=sce_atlas$sample)
  # genes_to_use <- rownames(decomp)[decomp$p.value<=0.01 & decomp$mean>0.1]
}

stopifnot(genes_to_use%in%rownames(sce_atlas))
stopifnot(genes_to_use%in%rownames(sce_query))

#########
## Map ##
#########

# TO-DO: TRY COSINE NORMALISATION
mapping  <- mapWrap(
  sce_atlas = sce_atlas,
  meta_atlas = meta_atlas,
  sce_query = sce_query,
  meta_query = meta_query,
  genes = genes_to_use,
  npcs = args$npcs,
  k = args$n_neighbours,
  cosineNorm = args$cosine_normalisation,
  order = NULL
)

################
## joint UMAP ##
################

# umap.mtx <- uwot::umap(rbind(mapping$pca_atlas_corrected, mapping$pca_query_corrected), n_neighbors=25, min_dist=0.50, metric="cosine")
# umap.dt <- umap.mtx %>% round(3) %>% as.data.table(keep.rownames = T) %>% setnames("rn","cell")
# 
# tmp <- rbind(
#   meta_atlas[,c("cell","celltype")] %>% .[,class:="atlas"],
#   mapping$mapping[,c("cell","celltype.mapped")] %>% as.data.table %>% setnames("celltype.mapped","celltype") %>% .[,class:="query"]
# )
# 
# to.plot <- umap.dt %>%
#   .[,cell:=c(meta_atlas$cell,meta_query$cell)] %>%
#   merge(tmp,by="cell")
# 
# to.plot.subset <- rbind(
#   to.plot[class=="query"],
#   to.plot[class=="atlas"][sample.int(.N, size=5e4)]
# ) %>% .[,class:=factor(class,levels=c("query","atlas"))] %>% setorder(-class)
# 
# p1 <- ggplot(to.plot.subset, aes(x=V1, y=V2, fill=class, size=class, alpha=class)) +
#   ggrastr::geom_point_rast(shape=21, stroke=0.1, raster.dpi=150) +
#   scale_size_manual(values=c("query"=1, "atlas"=0.5)) +
#   scale_alpha_manual(values=c("query"=0.9, "atlas"=0.65)) +
#   scale_fill_manual(values=c("atlas"="gray60", "query"="red")) +
#   theme_classic() +
#   ggplot_theme_NoAxes() +
#   theme(
#     legend.position="none"
#   )
# 
# p2 <- ggplot(to.plot.subset, aes(x=V1, y=V2, fill=celltype, size=class)) +
#   ggrastr::geom_point_rast(size=1, shape=21, stroke=0.1) +
#   scale_fill_manual(values=opts$celltype.colors) +
#   scale_size_manual(values=c("query"=1, "atlas"=0.5)) +
#   theme_classic() +
#   ggplot_theme_NoAxes() +
#   theme(
#     legend.position="none"
#   )
# 
# p <- cowplot::plot_grid(plotlist=list(p1,p2), nrow=1)
# 
# pdf(file.path(io$basedir,"results/mapping/pdf/fig/mapping_joint_umap.pdf"), width=12, height=5)
# print(p)
# dev.off()


##########
## Save ##
##########

mapping.dt <- mapping$mapping %>% 
  .[,c("cell","celltype.mapped","celltype.score","closest.cell")] %>% 
  as.data.table

# outfile <- sprintf("%s/mapping_mnn_%s.txt.gz",args$outdir,paste(args$query_samples,collapse="-"))
# fwrite(mapping.dt, outfile, sep="\t")
fwrite(mapping.dt, args$outfile, sep="\t")

##########
## TEST ##
##########

# previous_mapping.dt <- fread("/Users/argelagr/data/gastrulation_histones/results/rna/mapping/old/sample_metadata_after_mapping.txt.gz") %>%
#   .[,c("cell","celltype.mapped_mnn","celltype.mapped_seurat")] %>% setnames(c("cell","old_mapping_mnn","old_mapping_seurat"))
# foo <- merge(previous_mapping.dt, mapping.dt[,c("cell","celltype.mapped")], by="cell")
# 