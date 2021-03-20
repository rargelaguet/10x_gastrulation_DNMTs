library(SingleCellExperiment)
sce <- readRDS("/Users/ricard/data/10x_gastrulation_DNMTs/processed/all_batches/old/SingleCellExperiment.rds")
logcounts(sce)
assayNames(sce)
assays(sce)[["logcounts"]] <- NULL

sce <- sce[,sce$stage!="E12.5"]
print(object.size(sce), units="auto")

saveRDS(sce,"/Users/ricard/data/10x_gastrulation_DNMTs/processed/all_batches/SingleCellExperiment.rds")
