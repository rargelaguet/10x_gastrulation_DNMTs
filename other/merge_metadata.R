source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")

sample_metadata.1 <- fread("/Users/ricard/data/10x_gastrulation_DNMTs/processed/sixth_batch/metadata.txt.gz")
sample_metadata.2 <- fread("/Users/ricard/data/10x_gastrulation_DNMTs/sample_metadata.txt.gz")

colnames(sample_metadata.1)
colnames(sample_metadata.2)

cols <- c("cell", "batch", "barcode", "nFeature_RNA", "nCount_RNA", "percent.mt", "pass_QC", "class", "stage", "celltype.mapped")

sample_metadata <- rbind(sample_metadata.1[,cols,with=F], sample_metadata.2[,cols,with=F])

fwrite(sample_metadata, "/Users/ricard/data/10x_gastrulation_DNMTs/sample_metadata.txt.gz", sep="\t", quote=F, na="NA")
