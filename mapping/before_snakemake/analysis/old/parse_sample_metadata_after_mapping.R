source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")

files <- c(
	"/Users/ricard/data/10x_gastrulation_DNMTs/results/mapping/second_batch/mapping_mnn_DNMT3B_HET.txt.gz",
	"/Users/ricard/data/10x_gastrulation_DNMTs/results/mapping/second_batch/mapping_mnn_DNMT3B_KO.txt.gz",
	"/Users/ricard/data/10x_gastrulation_DNMTs/results/mapping/third_batch/mapping_mnn_Dnmt3aHet_Dnmt3bKO.txt.gz",
	"/Users/ricard/data/10x_gastrulation_DNMTs/results/mapping/third_batch/mapping_mnn_Dnmt3aKO_Dnmt3bWT.txt.gz",
	"/Users/ricard/data/10x_gastrulation_DNMTs/results/mapping/third_batch/mapping_mnn_Dnmt3aWT_Dnmt3bWT.txt.gz"
)

dt <- files %>% map(~ fread(.)[!is.na(celltype.mapped),c("barcode","batch","celltype.mapped","stage.mapped","closest.cell","celltype.score","cellstage.score")]) %>% rbindlist %>%
  .[,cell2:=paste(batch,barcode,sep="_")] %>%
  .[,c("barcode","batch"):=NULL]
  # .[cell2%in%sample_metadata$cell2]

sample_metadata <- fread("/Users/ricard/data/10x_gastrulation_DNMTs/processed/merged/sample_metadata.txt.gz") %>%
  .[,c("cell2", "cell", "barcode.x", "batch.x", "stage", "class", "target", "nCount_RNA", "nFeature_RNA", "percent.mt", "doublet_score", "pass_QC")] %>%
  setnames(c("barcode.x","batch.x"),c("barcode","batch")) %>%
  merge(dt,by="cell2",all.x=T)

dput(colnames(sample_metadata))

sample_metadata %>%
  .[,cell:=NULL] %>%
  setnames("cell2","cell")
fwrite(sample_metadata, io$metadata, sep="\t", col.names=T, row.names=F, na="NA", quote=F)
