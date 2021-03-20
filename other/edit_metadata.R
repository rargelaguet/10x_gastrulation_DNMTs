source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")

sample_metadata <- fread(io$metadata)

# Rename samples
foo <- sample_metadata[,c("batch","class")] %>% unique %>% .[,sample:=paste(class,1:.N,sep="_"),by="class"]

sample_metadata2 <- sample_metadata %>% merge(foo,by=c("batch","class")) %>%
  .[,c("cell","barcode","batch","sample","class","stage","nFeature_RNA","nCount_RNA", "percent.mt", "pass_QC", "celltype.mapped", "celltype.score","closest.cell")]

fwrite(sample_metadata2, io$metadata, sep="\t", quote=F, na="NA")
