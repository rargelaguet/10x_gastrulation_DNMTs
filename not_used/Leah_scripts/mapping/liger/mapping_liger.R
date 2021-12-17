#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(MouseGastrulationData))
suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(SeuratWrappers))
#suppressPackageStartupMessages(library(future)); plan("multiprocess", workers = 4)
options(future.globals.maxSize = 200000 * 1024^2)

###################
## Parse options ##
###################
 
option_list = list(
    make_option(c("-A", "--atlas.RDS"), type="character", default=NULL, 
                help="path to the batch effect corrected mouse gastrulation 10x atlas Seurat RDS", metavar="character"),
    make_option(c("-a", "--atlas.metadata"), type="character", default=NULL, 
                help="path to the batch effect corrected mouse gastrulation 10x atlas metdata txt", metavar="character"),
    make_option(c("-Q", "--query.seurat"), type="character", default=NULL, 
                help="path to the query to be mapped (a Seurat file)", metavar="character"),
    make_option(c("-q", "--query.metadata"), type="character", default=NULL, 
                help="path to the sample metadata file for the query", metavar="character"),
    make_option(c("-e", "--experiment"), type="character", default=NULL, 
                help="name of experiment", metavar="character"),
    make_option(c("-M", "--query.mapping.RDS"), type="character", default=NULL, 
                help="path to how the query previously mapped", metavar="character"),
    make_option(c("-m", "--query.mapping.meta"), type="character", default=NULL, 
                help="path to how the query previously mapped", metavar="character"),
    make_option(c("-b", "--query_batches"), type="character", default=NULL, 
                help="names of batches in query", metavar="character"),
    make_option(c("-S", "--settings"), type="character", default=NULL, 
                help="path to settings file", metavar="character"),
    make_option(c("-f", "--functions"), type="character", default=NULL, 
                help="path to file with functions used", metavar="character"),
    make_option(c("-O", "--output.RDS"), type="character", default=NULL, 
                help="path to output RDS file", metavar="character"),
    make_option(c("-o", "--output.metadata"), type="character", default=NULL, 
                help="path to output sample metadata file", metavar="character"),
    make_option(c("-l", "--lineage.list"), type="character", default=NULL, 
                help="which lineages to use", metavar="character"),
    make_option(c("-L", "--lineage.name"), type="character", default=NULL, 
                help="name of the lineage subset to use", metavar="character"),
    make_option(c("-s", "--atlas_stages"), type="character", default=NULL, 
                help="stages of the atlas to map to", metavar="character"),
    make_option(c("-t", "--testing"), action="store_true", default=NULL, 
                help="if given, data will be subset, so code can be tested")
); 
 
opt_parser = OptionParser(option_list=option_list);
opts = parse_args(opt_parser);

if (is.null(opts$query.seurat)){
    print_help(opt_parser)
    stop("The path to a query Seurat file must be supplied.n", call.=FALSE)
} else if (is.null(opts$query.metadata)){
    print_help(opt_parser)
    stop("The path to query sample metadata must be supplied.n", call.=FALSE)
} else if (is.null(opts$settings)) {
    print_help(opt_parser)
    stop("A settings file must be supplied.n", call.=FALSE)
} else if (is.null(opts$functions)) {
    print_help(opt_parser)
    stop("A functions file must be supplied.n", call.=FALSE)
} else if (is.null(opts$output.RDS)) {
    print_help(opt_parser)
    stop("A path to an RDS output file must be supplied.n", call.=FALSE)
} else if (is.null(opts$output.metadata)) {
    print_help(opt_parser)
    stop("A path to a metadata output file must be supplied.n", call.=FALSE)
} else if (((!is.null(opts$lineage.list)) || (!is.null(opts$lineage.name))) && (is.null(opts$query.mapping.RDS))) {
    print_help(opt_parser)
    stop("A path to a previous mapping output RDS file must be supplied.n", call.=FALSE)
} else if (((!is.null(opts$lineage.list)) || (!is.null(opts$lineage.name))) && (is.null(opts$query.mapping.meta))) {
    print_help(opt_parser)
    stop("A path to a previous mapping output metadata file must be supplied.n", call.=FALSE)
}

source(opts$settings)
source(opts$functions)

if ((!is.null(io$testing)) && (io$testing != TRUE)){
    stop("In the settings file, testing must be NULL or TRUE.n", call.=FALSE)
} else if ((!is.null(io$corrected.atlas.metadata)) && (!is.string(io$atlas.metadata))){
    stop("In the settings file, corrected.atlas.metadata must be NULL or a path.n", call.=FALSE)
} else if ((!is.null(io$corrected.atlas.RDS)) && (!is.string(io$atlas.RDS))){
    stop("In the settings file, corrected.atlas.RDS must be NULL or a path.n", call.=FALSE)
} else if ((!is.null(io$atlas_stages)) && (!is.string(io$atlas_stages))){
    stop("In the settings file, atlas_stages must be NULL or a path.n", call.=FALSE)
} else if ((!is.null(opts$lineage.name)) && (!opts$lineage.name %in% names(io$subset_celltypes))){
    stop("The lineage name must be the name of one of the settings file classes.n", call.=FALSE)
}

if (isTRUE(io$testing)) opts$testing = TRUE
if (is.string(io$atlas.RDS)) opts$atlas.RDS <- io$atlas.RDS
if (is.string(io$atlas.metadata)) opts$atlas.metadata <- io$atlas.metadata
if (is.string(io$atlas_stages)) opts$atlas_stages <- io$atlas_stages
if (is.null(opts$query_batches)) opts$query_batches = "all"
if (isTRUE(opts$testing)) message("Doing mapping in test mode (subsetting cells)...")

if (is.string(opts$query_batches)) opts$query_batches <- str_trim(strsplit(opts$query_batches, ",")[[1]])
if (is.string(opts$atlas_stages)) opts$atlas_stages <- str_trim(strsplit(opts$atlas_stages, ",")[[1]])
if (is.string(opts$lineage.list)) opts$lineage.list <- str_trim(strsplit(opts$lineage.list, ",")[[1]])
if (is.string(opts$lineage.name)) opts$lineage.list <- io$subset_celltypes[opts$lineage.name][[1]]

################
## Load atlas ##
################

message("Loading atlas...")

if ((is.null(opts$atlas.RDS)) || (is.null(opts$atlas.metadata))){
    
    message("Downloading atlas from bioconda package...")
    
    # Load SingleCellExperiment
    suppressMessages(sce_atlas <- EmbryoAtlasData(type="processed"))
    
    sce_atlas <- sce_atlas[,sce_atlas$doublet==F & sce_atlas$stripped==F]
    
    # Load metadata
    meta_atlas <- colData(sce_atlas) %>% as.data.table
    
} else {
    
    # Load metadata
    if (file_ext(opts$atlas.metadata) == "rds") {
        meta_atlas <- readRDS(paste0(opts$atlas.metadata)) %>% as.data.table
    } else {
        meta_atlas <- fread(opts$atlas.metadata, sep="\t", na="NA", quote=F) %>% as.data.table
    }
    
    meta_atlas <- meta_atlas[meta_atlas$doublet==FALSE,]
    meta_atlas <- meta_atlas[meta_atlas$stripped==FALSE,]
    
    # Load SingleCellExperiment
    sce_atlas  <- readRDS(paste0(opts$atlas.RDS))[,meta_atlas$cell]
}

# Filter stages
if (!is.null(opts$atlas_stages)) {
  message(sprintf("Subsetting atlas stages to %s",paste(opts$atlas_stages, collapse=", ")))
  meta_atlas <- meta_atlas[stage%in%opts$atlas_stages]
}

# Filter lineages
if (!is.null(opts$lineage.list)) {
  message(sprintf("Subsetting atlas to %s",paste(opts$lineage.list, collapse=", ")))
  meta_atlas <- meta_atlas[meta_atlas$celltype %in% opts$lineage.list]
}

# Filter if test mode
if (isTRUE(opts$testing)) {
    
    n <- round(10000/length(unique(meta_atlas$sample)))
    sub = c()
    for (s in unique(meta_atlas$sample)){
        sub <- append(sub, meta_atlas[meta_atlas$sample==s,]$cell[1:n])
    }
    
    meta_atlas <- meta_atlas[meta_atlas$cell %in% sub,]; rm(sub)
}


sce_atlas <- sce_atlas[,meta_atlas$cell]
sce_atlas <- multiBatchNorm(sce_atlas, batch=meta_atlas$stage, preserve.single=TRUE)
seurat_atlas <- as.Seurat(x=sce_atlas, project='ATLAS')
seurat_atlas@meta.data$map <- "ATLAS"
stopifnot(all(meta_atlas$cell == colnames(seurat_atlas)))


################
## Load query ##
################

message("Loading query...")

# Load metadata
meta_query <- fread(paste0(opts$query.metadata))

# Filter batches
if (any(opts$query_batches != "all")) {
  message(sprintf("Subsetting query batches to %s",paste(opts$query_batches, collapse=", ")))
  meta_query <- meta_query[class%in%opts$query_batches]
}

# Filter lineages
if (!is.null(opts$lineage.list)) {
  message(sprintf("Subsetting query to %s",paste(opts$lineage.list, collapse=", ")))
  meta_map <- fread(opts$query.mapping.meta)
  meta_map <- meta_map[meta_map$celltype.mapped %in% opts$lineage.list]
  meta_query <- merge(meta_query, meta_map, by=c(intersect(colnames(meta_map),colnames(meta_query))))
}


if ((!is.null(opts$lineage.list))) {
  
  if (dim(meta_map)[1] < 50) {
      message("Too few cells to remap. Instead subsetting previous mapping...")

      out <- readRDS(paste0(opts$query.mapping.RDS))
      out$closest.cells <- NULL
      out$celltypes.mapped <- NULL
      out$cellstages.mapped <- NULL
      out$correct_atlas <- out$correct_atlas[meta_atlas$cell,]
      out$correct_map <- out$correct_map[meta_map$cell,]
      out$integrated_seurat <- out$integrated_seurat[,c(meta_atlas$cell,meta_query$cell)]
      out$mapping <- meta_map

      # save mapping in .rds format
      #saveRDS(out, paste0(opts$output.RDS))

      # save mapping in .txt format
      #fwrite(meta_map, paste0(opts$output.metadata), sep="\t", na="NA", quote=F)

      opt <- options(show.error.messages = FALSE)
      on.exit(options(opt))
      stop()
  }

}

# Filter if test mode
if (isTRUE(opts$testing)) {
    
    n <- round(5000/length(opts$query_batches))
    sub = c()
    for (b in unique(meta_query$batch)){
        sub <- append(sub, meta_query[meta_query$batch==b,]$cell[1:n])
    }
    
    meta_query <- meta_query[meta_query$cell %in% sub,]; rm(sub)
}

# Load Seurat object
seurat_query <- readRDS(paste0(opts$query.seurat))[,meta_query$cell]
seurat_query@meta.data$map <- "QUERY"


stopifnot(all(meta_query$cell == colnames(seurat_query)))


#####################
## Intersect genes ## 
#####################

message("Filtering genes...")

# Remove genes with low expression in both query and atlas (not batch specific, as I'd otherwise likely filter out timepoint specific genes)
tmp1 <- rowSums(seurat_atlas@assays$RNA@counts); tmp1 <- tmp1[tmp1>=io$min.counts.per.gene]
tmp2 <- rowSums(seurat_query@assays$RNA@counts); tmp2 <- tmp2[tmp2>=io$min.counts.per.gene]

genes <- intersect(names(tmp1), names(tmp2))
seurat_query  <- seurat_query[genes,]
seurat_atlas <- seurat_atlas[genes,]


###########################
## Merge query and atlas ##
###########################

message("Merging query and atlas...")

seurat_merged <- merge(x=seurat_query, y=seurat_atlas, project="merged")
seurat_merged@meta.data$batch[which(is.na(seurat_merged@meta.data$batch))] <- seurat_merged@meta.data$sample[which(is.na(seurat_merged@meta.data$batch))] #atlas doesn't have a batch, called sample instead


#############################
## Normalise and find HVGs ##
#############################

message("Pre-processing data...")

# Basing code on: https://htmlpreview.github.io/?https://github.com/satijalab/seurat.wrappers/blob/master/docs/liger.html

seurat_merged <- NormalizeData(seurat_merged)
seurat_merged <- FindVariableFeatures(seurat_merged)
seurat_merged <- ScaleData(seurat_merged, split.by = "batch", do.center = FALSE)


##########################
## Integrating Datasets ##
##########################

message("Integrating datasets...")

seurat_merged <- RunOptimizeALS(seurat_merged, k = io$k, split.by = "batch")
seurat_merged <- RunQuantileAlignSNF(seurat_merged, split.by = "batch")


##############################
## Dimensionality Reduction ##
##############################

seurat_merged <- RunUMAP(seurat_merged, dims = 1:ncol(seurat_merged[["iNMF"]]), reduction = "iNMF", umap.method = "umap-learn")


############################
## Find Nearest Neighbors ##
############################

message("Find Nearest Neighbours...") 

md <- seurat_merged@meta.data
correct_map <- seurat_merged@reductions$iNMF@cell.embeddings[rownames(md)[which(md$map=='QUERY')],]
correct_atlas <- seurat_merged@reductions$iNMF@cell.embeddings[rownames(md)[which(md$map=='ATLAS')],]
atlas_meta <- md[which(md$map=='ATLAS'),]
map_meta <- md[which(md$map=='QUERY'),]
mapping <- get_meta(correct_atlas = correct_atlas,
                    atlas_meta = atlas_meta,
                    correct_map = correct_map,
                    map_meta = map_meta,
                    k_map = io$k)
message("Done\n")


#######################
## Get Mapping Score ##
#######################

message("Computing mapping scores...") 
out <- list()
for (i in seq(from = 1, to = io$k)) {
out$closest.cells[[i]]     <- sapply(mapping, function(x) x$cells.mapped[i])
out$celltypes.mapped[[i]]  <- sapply(mapping, function(x) x$celltypes.mapped[i])
out$cellstages.mapped[[i]] <- sapply(mapping, function(x) x$stages.mapped[i])
}  
multinomial.prob <- getMappingScore(out)
message("Done\n")

message("Writing output...") 
out$correct_atlas <- correct_atlas; rm(correct_atlas)
out$correct_map <- correct_map; rm(correct_map)
out$integrated_seurat <- seurat_merged; rm(seurat_merged)
ct <- sapply(mapping, function(x) x$celltype.mapped); is.na(ct) <- lengths(ct) == 0
st <- sapply(mapping, function(x) x$stage.mapped); is.na(st) <- lengths(st) == 0
cm <- sapply(mapping, function(x) x$cells.mapped[1]); is.na(cm) <- lengths(cm) == 0
out$mapping <- data.frame(
    cell            = names(mapping), 
    celltype.mapped = unlist(ct),
    stage.mapped    = unlist(st),
    closest.cell    = unlist(cm),
    stringsAsFactors=FALSE); rm(mapping)

out$mapping <- cbind(out$mapping,multinomial.prob)
message("Done\n")

##########
## Save ##
##########

# save mapping in .rds format
saveRDS(out, paste0(opts$output.RDS))

# save mapping in .txt format
metadata <- fread(paste0(opts$query.metadata)) %>%
  merge(out$mapping %>% as.data.table, by="cell", all.x=TRUE)
fwrite(metadata, paste0(opts$output.metadata), sep="\t", na="NA", quote=F)
