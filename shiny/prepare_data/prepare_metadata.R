source(here::here("settings.R"))

#####################
## Define settings ##
#####################

io$metadata <- file.path(io$basedir,"results/mapping/sample_metadata_after_mapping.txt.gz")
io$outfile <- paste0(io$basedir,"/shiny/cell_metadata.txt.gz")

opts$classes <- c("WT", "Dnmt3a_KO", "Dnmt3b_KO", "Dnmt1_KO")

opts$rename_celltypes <- c(
  "Erythroid3" = "Erythroid",
  "Erythroid2" = "Erythroid",
  "Erythroid1" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors",
  "Anterior_Primitive_Streak" = "Primitive_Streak"
)

########################
## Load cell metadata ##
########################

cell_metadata.dt <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE & class%in%opts$classes & !is.na(celltype.mapped)] %>%
  .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,opts$rename_celltypes)] %>%
  .[,dataset:=ifelse(grepl("Grosswendt",sample),"CRISPR","KO")] %>%
  .[,c("cell", "celltype.mapped","class", "alias","dataset","closest.cell")] %>%
  setnames(c("celltype.mapped","alias"),c("celltype","sample")) %>%
  .[,class_celltype:=sprintf("%s-%s",class,celltype)] %>%
  .[,class_celltype_dataset:=sprintf("%s-%s-%s",class,celltype,dataset)] %>%
  .[,class_sample_celltype_dataset:=sprintf("%s-%s-%s-%s",class,sample,celltype,dataset)]# %>%

table(cell_metadata.dt$sample)
table(cell_metadata.dt$class)
table(cell_metadata.dt$celltype)
table(cell_metadata.dt$dataset)

colnames(cell_metadata.dt)

stopifnot(sum(is.na(cell_metadata.dt$celltype))==0)

##########
## Save ##
##########

fwrite(cell_metadata.dt, io$outfile, sep="\t", na="NA", quote=F)
