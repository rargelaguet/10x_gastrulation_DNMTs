matrix.please<-function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[[1]]
  m
}

load_Seurat <- function(file, assay = "RNA", normalise = FALSE, features = NULL, cells = NULL, remove_non_expressed_genes = FALSE, ...) {
  library(Seurat)
  seurat <- readRDS(file)
  # if (assay%in%Seurat::Assays(seurat)) seurat <- seurat[[assay]]
  if (!is.null(cells)) seurat <- seurat[,cells]
  if (!is.null(features)) seurat <- seurat[features,]
  if (normalise) {
    seurat <- NormalizeData(seurat, normalization.method = "LogNormalize")
    seurat <- ScaleData(seurat, ...)
  }
  if (remove_non_expressed_genes) seurat <- seurat[which(Matrix::rowMeans(seurat@assays[[assay]]@counts)>1e-4),]
  return(seurat)
}

load_SingleCellExperiment <- function(file, normalise = FALSE, features = NULL, cells = NULL, remove_non_expressed_genes = FALSE) {
  library(SingleCellExperiment); library(scran); library(scater);
  sce <- readRDS(file)
  if (!is.null(cells)) {
    if (sum(!cells%in%colnames(sce)>0)) {
      message(sprintf("%s not found in the SingleCellExperiment object",paste(cells[!cells%in%colnames(sce)], collapse=", ")))
      cells <- cells[cells%in%colnames(sce)]
    }
    sce <- sce[,cells]
  }
  
  if (!is.null(features)) {
    if (sum(!features%in%rownames(sce)>0)) {
      message(sprintf("%s not found in the SingleCellExperiment object",paste(features[!features%in%rownames(sce)], collapse=", ")))
      features <- features[features%in%rownames(sce)]
    }
    sce <- sce[features,]
  }
  
  if (remove_non_expressed_genes) {
    sce <- sce[which(Matrix::rowSums(counts(sce))>15),]
  }
  
  if (normalise) {
    sce <- logNormCounts(sce)
  }
  
  return(sce)
}

linMap <- function(x, from, to) return( (x - min(x)) / max(x - min(x)) * (to - from) + from )

ggplot_theme_NoAxes <- function() {
  theme(
    axis.line = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
}

minmax.normalisation <- function(x)
{
  return((x-min(x,na.rm=T)) /(max(x,na.rm=T)-min(x,na.rm=T)))
}

# Remove unwanted effects from a matrix
#
# @parm mtx An expression matrix to regress the effects of covariates out
# of should be the complete expression matrix in genes x cells
# @param covariates A matrix or data.frame of latent variables, should be cells
# x covariates, the colnames should be the variables to regress
# @param features_idx An integer vector representing the indices of the
# genes to run regression on
# @param model.use Model to use, one of 'linear', 'poisson', or 'negbinom'; pass
# NULL to simply return mtx
# @param verbose Display a progress bar
#' @importFrom stats as.formula lm
#' @importFrom utils txtProgressBar setTxtProgressBar
#
RegressOutMatrix_parallel <- function(mtx, covariates = NULL, features_idx = NULL, split.by = NULL, block.size = 1000, min.cells.to.block = 3000, ncores = 1, verbose = TRUE) {
  
  library(future)
  library(future.apply)
  plan("multiprocess", workers = ncores)
  
  # Check features_idx
  if (is.null(features_idx)) {
    features_idx <- 1:nrow(mtx)
  }
  if (is.character(features_idx)) {
    features_idx <- intersect(features_idx, rownames(mtx))
    if (length(features_idx) == 0) {
      stop("Cannot use features that are beyond the scope of mtx")
    }
  } else if (max(features_idx) > nrow(mtx)) {
    stop("Cannot use features that are beyond the scope of mtx")
  }
  
  # Check dataset dimensions
  if (nrow(covariates) != ncol(mtx)) {
    stop("Uneven number of cells between latent data and expression data")
  }
  
  # Subset
  mtx <- mtx[features_idx,]
  mtx.dimnames <- dimnames(mtx)
  
  # Define chunck points
  chunk.points <- ChunkPoints(dsize = nrow(mtx), csize = block.size)
  
  # Define cell splitting
  split.cells <- split(colnames(mtx), f = split.by %||% TRUE)
  
  if (nbrOfWorkers() > 1) {
    
    # Define chuncks
    chunks <- expand.grid(
      names(split.cells),
      1:ncol(chunk.points),
      stringsAsFactors = FALSE
    )
    
    # Run RegressOutMatrix in parallel
    mtx.resid <- future_lapply(
      X = 1:nrow(chunks),
      FUN = function(i) {
        row <- chunks[i, ]
        group <- row[[1]]
        index <- as.numeric(row[[2]])
        return(RegressOutMatrix(
          mtx = mtx[chunk.points[1, index]:chunk.points[2, index], split.cells[[group]], drop = FALSE],
          covariates = covariates[split.cells[[group]], , drop = FALSE],
          # features_idx = features_idx[chunk.points[1, index]:chunk.points[2, index]],
          verbose = FALSE
        ))
      }
    )
    
    # Merge splitted cells
    if (length(split.cells) > 1) {
      merge.indices <- lapply(
        X = 1:length(x = split.cells),
        FUN = seq.int,
        to = length(mtx.resid),
        by = length(split.cells)
      )
      mtx.resid <- lapply(
        X = merge.indices,
        FUN = function(x) {
          return(do.call( 'rbind', mtx.resid[x]))
        }
      )
      mtx.resid <- do.call('cbind', mtx.resid)
    } else {
      mtx.resid <- do.call( 'rbind', mtx.resid)
    }
  } else {
    
    mtx.resid <- lapply(
      X = names(split.cells),
      FUN = function(x) {
        if (verbose && length(split.cells) > 1) {
          message("Regressing out variables from split ", x)
        }
        return(RegressOutMatrix(
          mtx = mtx[, split.cells[[x]], drop = FALSE],
          covariates = covariates[split.cells[[x]], , drop = FALSE],
          features_idx = features_idx,
          verbose = verbose
        ))
      }
    )
    mtx.resid <- do.call('cbind', mtx.resid)
  }
  # dimnames(mtx.resid) <- dimnames(mtx)
  return(mtx.resid)
}

RegressOutMatrix <- function(mtx, covariates = NULL, features_idx = NULL, verbose = TRUE) {
  
  # Check features_idx
  if (is.null(features_idx)) {
    features_idx <- 1:nrow(mtx)
  }
  if (is.character(features_idx)) {
    features_idx <- intersect(features_idx, rownames(mtx))
    if (length(features_idx) == 0) {
      stop("Cannot use features that are beyond the scope of mtx")
    }
  } else if (max(features_idx) > nrow(mtx)) {
    stop("Cannot use features that are beyond the scope of mtx")
  }
  
  # Check dataset dimensions
  if (nrow(covariates) != ncol(mtx)) {
    stop("Uneven number of cells between latent data and expression data")
  }
  
  # Subset
  mtx <- mtx[features_idx,]
  mtx.dimnames <- dimnames(mtx)
  
  # Create formula for regression
  vars.to.regress <- colnames(covariates)
  fmla <- paste('GENE ~', paste(vars.to.regress, collapse = '+')) %>% as.formula
  
  # In this code, we'll repeatedly regress different Y against the same X
  # (covariates) in order to calculate residuals.  Rather that repeatedly
  # call lm to do this, we'll avoid recalculating the QR decomposition for the
  # covariates matrix each time by reusing it after calculating it once
  regression.mat <- cbind(covariates, mtx[1,])
  colnames(regression.mat) <- c(colnames(covariates), "GENE")
  qr <- lm(fmla, data = regression.mat, qr = TRUE)$qr
  rm(regression.mat)
  
  # Make results matrix
  data.resid <- matrix(
    nrow = nrow(mtx),
    ncol = ncol(mtx)
  )
  
  if (verbose) pb <- txtProgressBar(char = '=', style = 3, file = stderr())
  
  # Extract residuals from each feature by using the pre-computed QR decomposition
  for (i in 1:length(features_idx)) {
    regression.mat <- cbind(covariates, mtx[features_idx[i], ])
    colnames(regression.mat) <- c(vars.to.regress, 'GENE')
    regression.mat <- qr.resid(qr = qr, y = mtx[features_idx[i],])  # The function qr.resid returns the residuals when fitting y to the matrix with QR decomposition.
    data.resid[i, ] <- regression.mat
    if (verbose) {
      setTxtProgressBar(pb = pb, value = i / length(features_idx))
    }
  }
  
  if (verbose) close(con = pb)
  
  dimnames(data.resid) <- mtx.dimnames
  
  return(data.resid)
}



# Generate chunk points
#
# @param dsize How big is the data being chunked
# @param csize How big should each chunk be
#
# @return A matrix where each column is a chunk, row 1 is start points, row 2 is end points
#
ChunkPoints <- function(dsize, csize) {
  return(vapply(
    X = 1L:ceiling(dsize / csize),
    FUN = function(i) {
      return(c(
        start = (csize * (i - 1L)) + 1L,
        end = min(csize * i, dsize)
      ))
    },
    FUN.VALUE = numeric(length = 2L)
  ))
}




smoother_aggregate_nearest_nb <- function(mat, D, k){
  # @param mat A matrix in a shape of #genes x #samples.
  # @param D A predefined distance matrix in a shape of #samples x #samples.
  # @param k An integer to choose \code{k} nearest samples (self-inclusive) to
  #  aggregate based on the distance matrix \code{D}. If \code{k} is greater than
  #  #samples, \code{k} is forced to be #samples to continue aggregation.
  sapply(seq_len(ncol(mat)), function(cid){
    nb_cid <- head(order(D[cid, ]), k)
    closest_mat <- mat[, nb_cid, drop=FALSE]
    # return(Matrix::rowSums(closest_mat))
    return(Matrix::rowMeans(closest_mat))
  })
}

give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}

sort.abs <- function(dt, sort.field) dt[order(-abs(dt[[sort.field]]))]


#' Write count data in the 10x format
#' Create a directory containing the count matrix and cell/gene annotation from a sparse matrix of UMI counts, 
#' in the format produced by the CellRanger software suite.
# https://rdrr.io/github/MarioniLab/DropletUtils/src/R/write10xCounts.R
write10xCounts <- function(path, x, barcodes=colnames(x), gene.id=rownames(x), gene.symbol=gene.id, gene.type="Gene Expression",
                           overwrite=FALSE, version=c("2", "3"))
{
  # Doing all the work on a temporary location next to 'path', as we have permissions there.
  # This avoids problems with 'path' already existing.
  temp.path <- tempfile(tmpdir=dirname(path)) 
  on.exit({ 
    if (file.exists(temp.path)) { unlink(temp.path, recursive=TRUE) } 
  })
  
  # Checking the values.
  if (length(gene.id)!=length(gene.symbol) || length(gene.id)!=nrow(x)) {
    stop("lengths of 'gene.id' and 'gene.symbol' must be equal to 'nrow(x)'")
  }
  if (ncol(x)!=length(barcodes)) { 
    stop("'barcodes' must of of the same length as 'ncol(x)'")
  }
  
  # Determining what format to save in.
  version <- match.arg(version)
  .write_sparse(temp.path, x, barcodes, gene.id, gene.symbol, gene.type, version=version)
  
  # We don't put this at the top as the write functions might fail; 
  # in which case, we would have deleted the existing 'path' for nothing.
  if (overwrite) {
    unlink(path, recursive=TRUE)
  } else if (file.exists(path)) { 
    stop("specified 'path' already exists")
  }
  file.rename(temp.path, path)
  return(invisible(TRUE))
}

.write_sparse <- function(path, x, barcodes, gene.id, gene.symbol, gene.type, version="2") {
  dir.create(path, showWarnings=FALSE)
  # gene.info <- data.frame(gene.id, gene.symbol, stringsAsFactors=FALSE)
  gene.info <- data.frame(gene.symbol, gene.id, stringsAsFactors=FALSE)
  
  if (version=="3") {
    gene.info$gene.type <- rep(gene.type, length.out=nrow(gene.info))
    mhandle <- file.path(path, "matrix.mtx")
    bhandle <- gzfile(file.path(path, "barcodes.tsv.gz"), open="wb")
    fhandle <- gzfile(file.path(path, "features.tsv.gz"), open="wb")
    on.exit({
      close(bhandle)
      close(fhandle)
    })
  } else {
    mhandle <- file.path(path, "matrix.mtx")
    bhandle <- file.path(path, "barcodes.tsv")
    fhandle <- file.path(path, "genes.tsv")
  }
  
  Matrix::writeMM(x, file=mhandle)
  write(barcodes, file=bhandle)
  utils::write.table(gene.info, file=fhandle, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  
  if (version=="3") {
    # Annoyingly, writeMM doesn't take connection objects.
    R.utils::gzip(mhandle)
  }
  
  return(NULL)
}


