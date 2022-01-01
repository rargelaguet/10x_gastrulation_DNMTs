# function to pseudobulk a SingleCellExperiment object
pseudobulk_sce_fn <- function(x, assay = NULL, by, fun = c("sum", "mean", "median"), scale = FALSE) {
  
  # check validity of input arguments
  fun <- match.arg(fun)
  if (is.null(assay))  assay <- assayNames(x)[1] 
  
  # store aggregation parameters &
  # nb. of cells that went into aggregation
  md <- metadata(x)
  md$agg_pars <- list(assay = assay, by = by, fun = fun, scale = scale)
  
  # get aggregation function
  # fun <- switch(fun, sum = "rowSums", mean = "rowMeans", median = "rowMedians")
  
  # drop missing factor levels & tabulate number of cells
  cd <- dplyr::mutate_if(as.data.frame(colData(x)), is.factor, droplevels)
  colData(x) <- DataFrame(cd, row.names = colnames(x),check.names = FALSE)
  md$n_cells <- table(as.data.frame(colData(x)[, by]))
  
  # assure 'by' colData columns are factors so that missing combinations aren't dropped
  for (i in by) 
    if (!is.factor(x[[i]])) 
      x[[i]] <- factor(x[[i]])
  
  # split cells & compute pseudo-bulks
  cs <- .split_cells(x, by)
  # pb <- .pb(x, cs, assay, fun)
  pb <- .pb(x=x, by=by, fun=fun)
  if (scale & length(by) == 2) {
    ls <- lapply(.pb(x, cs, "counts", "rowSums"), colSums)
    pb <- lapply(seq_along(pb), function(i) pb[[i]] / 1e6 * ls[[i]])
    names(pb) <- names(ls)
  }
  
  # construct SCE
  pb <- SingleCellExperiment(pb, metadata = md)
  
  # propagate 'colData' columns that are unique across 2nd 'by'
  if (length(by) == 2) {
    cd <- colData(x)
    ids <- colnames(pb)
    counts <- vapply(ids, function(u) {
      m <- as.logical(match(cd[, by[2]], u, nomatch = 0))
      vapply(cd[m, ], function(u) length(unique(u)), numeric(1))
    }, numeric(ncol(colData(x))))
    cd_keep <- apply(counts, 1, function(u) all(u == 1))
    cd_keep <- setdiff(names(which(cd_keep)), by)
    if (length(cd_keep) != 0) {
      m <- match(ids, cd[, by[2]], nomatch = 0)
      cd <- cd[m, cd_keep, drop = FALSE]
      rownames(cd) <- ids
      colData(pb) <- cd
    }
  }
  return(pb)
}


# split cells by cluster-sample
# auxiliary function to pseudobulk a SingleCellExperiment object
#   - by: character vector specifying colData column(s) to split by. If length(by) == 1, a list of length nlevels(colData$by), else,
#          a nested list with 2nd level of length nlevels(colData$by[2])
.split_cells <- function(x, by) {
  if (is(x, "SingleCellExperiment")) x <- colData(x)
  cd <- data.frame(x[by], check.names = FALSE)
  cd <- data.table(cd, cell = rownames(x)) %>% split(by = by, sorted = TRUE, flatten = FALSE)
  purrr::map_depth(cd, length(by), "cell")
}


# auxiliary function to pseudobulk a SingleCellExperiment object
.pb <- function(x, by, fun) {
  
  # compute pseudobulks
  # y <- scuttle::summarizeAssayByGroup(x, assay.type = assay, ids = (ids <- colData(x)[by]), statistics = fun, BPPARAM = BiocParallel::SerialParam())
  y <- scuttle::summarizeAssayByGroup(x, ids = colData(x)[by], statistics = fun)
  colnames(y) <- y[[by[length(by)]]]
  
  if (length(by) == 1)  return(assay(y))
  
  # reformat into one assay per 'by[1]'
  if (is.factor(ids <- y[[by[1]]]))
    ids <- droplevels(ids)
  is <- split(seq_len(ncol(y)), ids)
  ys <- map(is, ~assay(y)[, .])
  
  # fill in missing combinations
  for (i in seq_along(ys)) {
    fill <- setdiff(unique(y[[by[2]]]), colnames(ys[[i]]))
    if (length(fill != 0)) {
      foo <- matrix(0, nrow(x), length(fill))
      colnames(foo) <- fill
      foo <- cbind(ys[[i]], foo)
      o <- paste(sort(unique(y[[by[2]]])))
      ys[[i]] <- foo[, o]
    }
  }
  return(ys)
}