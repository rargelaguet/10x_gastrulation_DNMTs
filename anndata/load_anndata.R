##########
## TEST ##
##########

library(zellkonverter)
# h5ad <- system.file("extdata", "krumsiek11.h5ad", package = "zellkonverter")
h5ad <- "/Users/ricard/data/gastrulation10x/processed/scanpy/Embryo10Xv6_genes_PCAbatchCorrected_mergedAGA.h5ad"
sce <- readH5AD(h5ad, use_hdf5 = T)



X.delayed <- sce@assays@data$X
print(object.size(X.delayed), units="auto")

X.matrix <- as.matrix(sce@assays@data$X)
print(object.size(X.matrix), units="auto")

print(object.size(sce), units="auto")
