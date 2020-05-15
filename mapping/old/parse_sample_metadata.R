
######################
## Define  settings ##
######################

source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")

# I/O
io$seurat <- "/Users/ricard/data/10x_gastrulation_DNMTs/processed/seurat.rds"
io$mapping <- "/Users/ricard/data/10x_gastrulation_DNMTs/mapping/mapping10x_mnn.rds"
io$outdir <- "/Users/ricard/data/10x_gastrulation_DNMTs/mapping"

###############
## Load data ##
###############

# seurat <- readRDS(io$seurat)
mapping <- readRDS(io$mapping)

sample_metadata <- seurat@meta.data %>% as.data.table %>%
  .[cell%in%mapping$mapping$cell]

################################################
## Add lineage annotations to sample metadata ##
################################################

mapping.dt <- mapping$mapping %>% as.data.table %>% 
  setnames("celltype.mapped", "celltype.mapped.level1") %>%
  .[,closest.cell:=NULL]

sample_metadata <- merge(sample_metadata,mapping.dt, by="cell", all.x=T)

sample_metadata %>%
  .[,celltype.mapped.level2:=celltype.mapped.level1] %>%
  
  # Mesoderm
  .[celltype.mapped.level1%in%c("Pharyngeal mesoderm","Paraxial mesoderm","Mixed mesoderm","Intermediate mesoderm","Somitic mesoderm","Caudal mesoderm", "Caudal Mesoderm"), celltype.mapped.level2:="Mesoderm"] %>%
  # Blood
  .[celltype.mapped.level1%in%c("Erythroid1","Erythroid2","Erythroid3"), celltype.mapped.level2:="Erythroid"] %>%
  .[celltype.mapped.level1%in%c("Blood progenitors 1","Blood progenitors 2"), celltype.mapped.level2:="Blood progenitors"] %>%
  # Embryonic endoderm
  .[celltype.mapped.level1%in%c("Gut","Def. endoderm"), celltype.mapped.level2:="Endoderm"] %>%
  # Extra-embryonic endoderm
  .[celltype.mapped.level1%in%c("Parietal endoderm","ExE endoderm","Visceral endoderm"), celltype.mapped.level2:="ExE endoderm"] %>%
  # Primitive streak
  .[celltype.mapped.level1%in%c("Primitive Streak", "Caudal epiblast","Anterior Primitive Streak"), celltype.mapped.level2:="Primitive Streak"] %>%
  # Ectoderm
  .[celltype.mapped.level1%in%c("Rostral neurectoderm","Surface ectoderm","Caudal neurectoderm"), celltype.mapped.level2:="Ectoderm"]

unique(sample_metadata$celltype.mapped.level2)

# Save
fwrite(sample_metadata, file=paste0(io$outdir,"/sample_metadata_mapping_mnn.txt"), sep="\t", row.names=F, col.names=T, na="NA", quote=F)

foo <- sample_metadata %>% as.data.frame %>% tibble::column_to_rownames("cell")
foo$cell <- rownames(foo)

# seurat@meta.data <- foo
# saveRDS(seurat, "/Users/ricard/data/10x_gastrulation_DNMTs/processed/seurat_2.rds")
