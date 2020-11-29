#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Seurat))
# suppressPackageStartupMessages(library(assertthat))
# suppressPackageStartupMessages(library(scran))


#####################
## Define settings ##
#####################

source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")

# Define options
opts$test <- TRUE

# Define I/O
io$inputdir <- paste0(io$basedir,"/original/all_batches")
io$outputdir <- paste0(io$basedir,"/processed/all_batches")

##############################
## Load and merge data sets ##
##############################

if (opts$test) opts$batches <- head(opts$batches,n=5)

mtx <- list()
cell.info <- list()
gene.info <- list()

for (i in opts$batches) {
  print(i)
    
  # Load gene metadata
  gene.loc <- sprintf("%s/%s_features.tsv.gz",io$inputdir,i)
  gene.info[[i]] <- read.delim(gene.loc, header = FALSE, colClasses = "character", stringsAsFactors = FALSE, sep="\t")[,c(1,2)]
  colnames(gene.info[[i]]) <- c("ens_id","symbol")
  rownames(gene.info[[i]]) <- NULL
  
  # Load cell metadata
  barcode.loc <- sprintf("%s/%s_barcodes.tsv.gz",io$inputdir,i)
  cell.info[[i]] <- read.delim(barcode.loc, header = FALSE, colClasses = "character", stringsAsFactors = FALSE, sep="\t")
  colnames(cell.info[[i]]) <- c("barcode")
  cell.info[[i]]$batch <- i
  cell.info[[i]]$cell <- sprintf("%s_%s",i,cell.info[[i]]$barcode)
  
  # Load matrix  
  matrix.loc <- sprintf("%s/%s_matrix.mtx.gz",io$inputdir,i)
  # matrix.loc <- sprintf("%s/%s/soup/soupX_adjusted_matrix.mtx.gz",io$inputdir,i)
  mtx[[i]] <- Matrix::readMM(matrix.loc)
  rownames(mtx[[i]]) <- gene.info[[i]]$symbol
  colnames(mtx[[i]]) <- cell.info[[i]]$barcode
}


# remove human genes for batches that are mapped to the mixed transcriptome
for (i in opts$batches) {
  if (any(grep("mm10",gene.info[[i]]$ens_id))) {
    idx <- grep("mm10",gene.info[[i]]$ens_id)
    gene.info[[i]] <- gene.info[[i]][idx,]
    gene.info[[i]]$ens_id <- gene.info[[i]]$ens_id %>% gsub("mm10___","",.)
    mtx[[i]] <- mtx[[i]][idx,]
    rownames(mtx[[i]]) <- rownames(mtx[[i]]) %>% gsub("mm10___","",.)
  }
}

# sanity checks
stopifnot(length(unique(lapply(mtx,nrow)))==1)
stopifnot(length(unique(lapply(mtx,rownames)))==1)

# bind gene names and remove human alignments
# gene.info <- do.call("rbind", gene.info)
# rownames(gene.info) <- NULL
# gene.info <- unique(gene.info)
# gene.info <- gene.info[grepl('mm10', gene.info$symbol),]
# gene.info$ens_id <- stringr::str_split_fixed(gene.info$ens_id,"___",2)[,2]
# gene.info$symbol <- stringr::str_split_fixed(gene.info$symbol,"___",2)[,2]
# rownames(gene.info) <- NULL

# Concatenate cell  metadata
cell.info <- do.call("rbind",cell.info)
rownames(cell.info) <- cell.info$cell

# Concatenate matrices
mtx <- do.call("cbind",mtx)
colnames(mtx) <- cell.info$cell
# mtx <- mtx[grepl('mm10', rownames(mtx)),]
# rownames(mtx) <- stringr::str_split_fixed(rownames(mtx),"___",2)[,2]

# Sanity checks
# cell.info$cell[!cell.info$cell %in% sample_metadata$cell]

################
## Processing ##
################

# Optionally subset protein-coding genes
# if (!is.null(opts$subset.proteincoding)){
#     genes <- fread(opts$subset.proteincoding)[,ens_id]
#     genes <- genes[genes %in% mouse.genes]
#     mouse.genes <- mouse.genes[mouse.genes %in% genes]
#     mtx <- mtx[mouse.genes,]
# }

# Subset cell metadata
# cell.info <- cell.info[colnames(mtx),]

# Subset gene metadata
gene.info <- gene.info[!duplicated(gene.info$symbol),]

# Subset matrix
mtx <- mtx[gene.info$symbol,]

# Sanity checks
sum(duplicated(rownames(mtx)))
sum(duplicated(colnames(mtx)))
stopifnot(all(colnames(mtx) == cell.info$cell))
stopifnot(all(rownames(mtx) == gene.info$symbol))

##########################
## Create Seurat object ##
##########################

srat <- CreateSeuratObject(mtx, meta.data = cell.info)

srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^mt-")

##########
## Save ##
##########

metadata <- srat@meta.data %>% as.data.table %>% .[,orig.ident:=NULL] %>%
  .[,c("cell","barcode","batch","nFeature_RNA","nCount_RNA","percent.mt")]

saveRDS(srat, paste0(io$outputdir,"/seurat.rds"))
fwrite(cell.info, paste0(io$outputdir,"/cell_info.txt.gz"), quote=F, na="NA", sep="\t")
fwrite(gene.info, paste0(io$outputdir,"/gene_info.txt.gz"), quote=F, na="NA", sep="\t")
fwrite(metadata, paste0(io$outputdir,"/metadata.txt.gz"), quote=F, na="NA", sep="\t")

