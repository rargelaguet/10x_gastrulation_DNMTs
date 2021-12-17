#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(scran))

###################
## Parse options ##
###################
 
option_list = list(
    make_option(c("-i", "--inputdir"), type="character", default=NULL, 
                help="directory of the input files", metavar="character"),
    make_option(c("-S", "--settings"), type="character", default=NULL, 
                help="path to settings file", metavar="character"),
    make_option(c("-O", "--output.seurat"), type="character", default=NULL, 
                help="path to output seurat file", metavar="character"),
    make_option(c("-o", "--output.metadata"), type="character", default=NULL, 
                help="path to output sample metadata file", metavar="character"),
    make_option(c("-g", "--output.gene.metadata"), type="character", default=NULL, 
                help="path to output gene metadata file", metavar="character"),
    make_option(c("-e", "--experiment"), type="character", default=NULL, 
                help="which experiment the data comes from", metavar="character"),
    make_option(c("-p", "--subset.proteincoding"), type="character", default=NULL, 
                help="path to proteincoding genes, if subsetting", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opts = parse_args(opt_parser);

if (is.null(opts$inputdir)){
    print_help(opt_parser)
    stop("An input directory must be supplied.n", call.=FALSE)
} else if (is.null(opts$settings)) {
    print_help(opt_parser)
    stop("A settings file must be supplied.n", call.=FALSE)
} else if (is.null(opts$output.seurat)) {
    print_help(opt_parser)
    stop("A directory to a Seurat output file must be supplied.n", call.=FALSE)
} else if (is.null(opts$output.metadata)) {
    print_help(opt_parser)
    stop("A directory to a sample metadata output file must be supplied.n", call.=FALSE)
} else if (is.null(opts$experiment)) {
    print_help(opt_parser)
    stop("The experiment must be specified.n", call.=FALSE)
}

##########################
## Source settings file ##
##########################

source(opts$settings)

if ((!is.null(io$qc)) && (io$qc != TRUE)){
    stop("In the settings file, qc must be NULL or TRUE.n", call.=FALSE)
} else if ((!is.null(io$subset.proteincoding)) && (!is.string(io$subset.proteincoding))){
    stop("In the settings file, subset.proteincoding must be NULL or a path.n", call.=FALSE)
} else if ((!is.null(io$scaledata)) && (io$scaledata != TRUE)){
    stop("In the settings file, scaledatav must be NULL or TRUE.n", call.=FALSE)
} else if ((opts$experiment != "second_batch") && (opts$experiment != "third_batch") && (opts$experiment != "fourth_batch")) {
    stop("The experiment must be either second_batch, third_batch, or fourth_batch.n", call.=FALSE)
}


##############################
## Load and merge data sets ##
##############################

message("Loading datasets")

mtx <- list()
cell.info <- list()
gene.info <- list()

toiter <- io$samples[[opts$experiment]]
for (i in toiter) {
    
  # Load cell metadata
  barcode.loc <- sprintf("%s%s_barcodes.tsv.gz",opts$inputdir,i)
  cell.info[[i]] <- read.delim(barcode.loc, header = FALSE, colClasses = "character", stringsAsFactors = FALSE, sep="\t")
  colnames(cell.info[[i]]) <- c("barcode")
  cell.info[[i]]$batch <- i
  tmp <- strsplit(i,"_")[[1]]
  
  if (opts$experiment == "second_batch") {
      
      cell.info[[i]]$stage <- tmp[1]
      if (grepl("Rep", tmp[2], fixed=TRUE)) {
          cell.info[[i]]$Rep <- tmp[2]
      } else {
          cell.info[[i]]$Rep <- NA
      }
      if (grepl("WT", i, fixed=TRUE)) {
          cell.info[[i]]$class <- "WT"
          cell.info[[i]]$target <- NA
          cell.info[[i]]$modification <- NA
      } else if (grepl("TET", i, fixed=TRUE)){
          cell.info[[i]]$class <- paste("TET",tmp[which(tmp=="TET")+1],sep="_")
          cell.info[[i]]$target <- "TET"
          cell.info[[i]]$modification <- tmp[which(tmp=="TET")+1]
      } else if (grepl("DNMT3A", i, fixed=TRUE)){
          # NOTE: RENAMING TO DNMT3B BECAUSE IT WAS A MISNAMING WHEN HANDLING THE FILES
          cell.info[[i]]$class <- paste("DNMT3B",tmp[which(tmp=="DNMT3A")+1],sep="_")
          cell.info[[i]]$target <- "DNMT3B"
          cell.info[[i]]$modification <- tmp[which(tmp=="DNMT3A")+1]
      } else {
          stop("Didn't recognise sample name based on hard-coding")
      }
      cell.info[[i]]$lane <- tmp[4]
      
  } else if (opts$experiment == "third_batch") {
      
    cell.info[[i]]$index <- tmp[1]
    cell.info[[i]]$stage <- tmp[2]
    cell.info[[i]]$embryo <- paste(tmp[2],tmp[3],sep="_")
    if (grepl("TET", i, fixed=TRUE)){
        cell.info[[i]]$class <- paste("TET",tmp[6],sep="_")
        cell.info[[i]]$target <- tmp[4]
        cell.info[[i]]$modification <- tmp[6]
        cell.info[[i]]$type <- tmp[5]
        cell.info[[i]]$Dnmt3a <- NA
        cell.info[[i]]$Dnmt3b <- NA
    } else if (grepl("Dnmt", i, fixed=TRUE)) {
        a_mod <- strsplit(tmp[4], "Dnmt3a")[[1]][2]
        cell.info[[i]]$class <- paste0(tmp[4],"_",tmp[5],tmp[6])
        cell.info[[i]]$target <- "Dnmt3a_Dnmt3b"
        cell.info[[i]]$modification <- paste0(a_mod,"_",tmp[6])
        cell.info[[i]]$type <- "DoubleMutant"
        cell.info[[i]]$Dnmt3a <- a_mod
        cell.info[[i]]$Dnmt3b <- tmp[6]
    }
    cell.info[[i]]$lane <- tmp[7]
  
  } else if (opts$experiment == "fourth_batch") {
  
    print("fourth_batch")
    tmp <- strsplit(toupper(i),"_")[[1]]
    cell.info[[i]]$lib_num <- tmp[1]
    cell.info[[i]]$stage <- paste(tmp[2],tmp[3],sep=".")
    cell.info[[i]]$target <- "Dnmt3a_Dnmt3b"
    if ((grepl("D3A", toupper(i), fixed=TRUE)) && (grepl("D3B", toupper(i), fixed=TRUE))) {
        cell.info[[i]]$class <- paste0(tmp[2],tmp[3],"_","Dnmt3a",tmp[which(tmp=="D3A")+1],"_","Dnmt3b",tmp[which(tmp=="D3B")+1])
        cell.info[[i]]$modification <- paste0(tmp[which(tmp=="D3A")+1],"_",tmp[which(tmp=="D3B")+1])
        if ((tmp[which(tmp=="D3A")+1] == "WT") && (tmp[which(tmp=="D3B")+1] == "WT")) {
            cell.info[[i]]$modification <- "WT"
        } else if ((tmp[which(tmp=="D3A")+1] == "WT") && (tmp[which(tmp=="D3B")+1] != "WT")) {
            cell.info[[i]]$modification <- "SingleMutant"
        } else if ((tmp[which(tmp=="D3A")+1] != "WT") && (tmp[which(tmp=="D3B")+1] == "WT")) {
            cell.info[[i]]$modification <- "SingleMutant"
        } else if ((tmp[which(tmp=="D3A")+1] != "WT") && (tmp[which(tmp=="D3B")+1] != "WT")) {
            cell.info[[i]]$modification <- "DoubleMutant"
        }
        cell.info[[i]]$Dnmt3a <- tmp[which(tmp=="D3A")+1]
        cell.info[[i]]$Dnmt3b <- tmp[which(tmp=="D3B")+1]
    } else if (grepl("D3A", toupper(i), fixed=TRUE)) {
        cell.info[[i]]$class <- paste0(tmp[2],tmp[3],"_","Dnmt3a",tmp[which(tmp=="D3A")+1],"_","Dnmt3b","WT")
        cell.info[[i]]$modification <- paste0(tmp[which(tmp=="D3A")+1],"_","WT")
        if (tmp[which(tmp=="D3A")+1] == "WT") {
            cell.info[[i]]$modification <- "WT"
        } else {
            cell.info[[i]]$modification <- "SingleMutant"
        }
        cell.info[[i]]$Dnmt3a <- tmp[which(tmp=="D3A")+1]
        cell.info[[i]]$Dnmt3b <- "WT"
    } else if (grepl("D3B", toupper(i), fixed=TRUE)) {
        cell.info[[i]]$class <- paste0(tmp[2],tmp[3],"_","Dnmt3a","WT","_","Dnmt3b",tmp[which(tmp=="D3B")+1])
        cell.info[[i]]$modification <- paste0("WT","_",tmp[which(tmp=="D3b")+1])
        if (tmp[which(tmp=="D3B")+1] == "WT") {
            cell.info[[i]]$modification <- "WT"
        } else {
            cell.info[[i]]$modification <- "SingleMutant"
        }
        cell.info[[i]]$Dnmt3a <- "WT"
        cell.info[[i]]$Dnmt3b <- tmp[which(tmp=="D3B")+1]
    }
  
  }
  
  # Load gene metadata (note we could just load this once)
  gene.loc <- sprintf("%s%s_features.tsv.gz",opts$inputdir,i)
  gene.info[[i]] <- read.delim(gene.loc, header = FALSE, colClasses = "character", stringsAsFactors = FALSE, sep="\t")[,c(1,2)]
  colnames(gene.info[[i]]) <- c("ens_id","symbol")
  
  # Load matrix  
  matrix.loc <- sprintf("%s%s_matrix.mtx.gz",opts$inputdir,i)
  mtx[[i]] <- Matrix::readMM(matrix.loc)
  rownames(mtx[[i]]) <- gene.info[[i]]$symbol
  colnames(mtx[[i]]) <- cell.info[[i]]$barcode
}

# bind gene names and remove human alignments
gene.info <- do.call("rbind", gene.info)
rownames(gene.info) <- NULL
gene.info <- unique(gene.info)
if (opts$experiment != "fourth_batch") {
    gene.info <- gene.info[grepl('mm10', gene.info$symbol),]
    gene.info$ens_id <- stringr::str_split_fixed(gene.info$ens_id,"___",2)[,2]
    gene.info$symbol <- stringr::str_split_fixed(gene.info$symbol,"___",2)[,2]
    rownames(gene.info) <- NULL
}

# Concatenate cell  metadata
cell.info <- do.call("rbind",cell.info)
cell.info$cell <- paste("cell",1:nrow(cell.info),sep="_")
rownames(cell.info) <- cell.info$cell

# Concatenate matrices
mtx <- do.call("cbind",mtx)
colnames(mtx) <- cell.info$cell
if (opts$experiment != "fourth_batch") {
    mtx <- mtx[grepl('mm10', rownames(mtx)),]
    rownames(mtx) <- stringr::str_split_fixed(rownames(mtx),"___",2)[,2]
}

################
## Processing ##
################

message("Processing datasets...")

# Optionally subset protein-coding genes
if (!is.null(opts$subset.proteincoding)){
    genes <- fread(opts$subset.proteincoding)[,ens_id]
    genes <- genes[genes %in% mouse.genes]
    mouse.genes <- mouse.genes[mouse.genes %in% genes]
    mtx <- mtx[mouse.genes,]
}

# Subset cell metadata
cell.info <- cell.info[colnames(mtx),]

########
## QC ##
########

if (io$qc == TRUE) {

    message("QCing...")

    srat <- CreateSeuratObject(mtx, meta.data = cell.info)

    srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^mt-")
    md <- srat@meta.data
    md$doublet_score <- NA

    mtx_sub_list <- list()
    cell.info_sub_list <- list()
    srat_sub_list <- list()
    batches <- unique(cell.info$batch)
    for (b in batches) {
        foo <- io$min_nFeature_RNA[[opts$experiment]][[b]]
        bar <- io$min_nCount_RNA[[opts$experiment]][[b]]
        baz <- io$max_percent.mt[[opts$experiment]][[b]]
        srat_sub_list[[b]] <- subset(srat, subset = nFeature_RNA > foo & nCount_RNA > bar & percent.mt < baz, cells = rownames(srat@meta.data)[which(srat@meta.data$batch == b)])
        sce_tmp <- as.SingleCellExperiment(srat_sub_list[[b]])
        srat_sub_list[[b]]@meta.data$doublet_score <- doubletCells(sce_tmp)
        md[rownames(srat_sub_list[[b]]@meta.data),"doublet_score"] <- srat_sub_list[[b]]@meta.data$doublet_score
        srat_sub_list[[b]] <- subset(srat_sub_list[[b]], subset = doublet_score < io$max_doublet_score[[opts$experiment]][[b]])
        cell.info_sub_list[[b]] <- srat_sub_list[[b]]@meta.data
        cell.info_sub_list[[b]]$origrn <- rownames(cell.info_sub_list[[b]])
        mtx_sub_list[[b]] <- srat_sub_list[[b]]@assays$RNA@counts
    }
    mtx <- do.call("cbind", mtx_sub_list)
    old_rownames <- rownames(mtx)
    rownames(mtx) <- gene.info$ens_id
    new_rownames <- rownames(mtx)
    tmp.cell.info <- do.call("rbind", cell.info_sub_list)
    rownames(tmp.cell.info) <- tmp.cell.info$origrn
    md$pass_QC <- FALSE
    md[rownames(tmp.cell.info),"pass_QC"] <- TRUE
    cell.info <- md
    rm(sce_tmp, srat_sub_list, cell.info_sub_list, mtx_sub_list, srat, md)
}

############
## Seurat ##
############

# Create seurat object

message("Creating seurat...")

srat <- CreateSeuratObject(mtx, meta.data = cell.info)


##########################
## SingleCellExperiment ##
##########################

# sce <- as.SingleCellExperiment(srat)

##########
## Save ##
##########

saveRDS(srat, opts$output.seurat)
fwrite(cell.info, opts$output.metadata, quote=F, na="NA", sep="\t")
fwrite(gene.info, opts$output.gene.metadata, quote=F, na="NA", sep="\t")
