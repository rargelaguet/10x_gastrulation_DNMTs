source(here::here("settings.R"))
#####################
## Define settings ##
#####################

io$metadata <- file.path(io$basedir,"results/mapping/sample_metadata_after_mapping.txt.gz")
io$outfile <- "/Users/argelagr/data/shiny_dnmt_tet/cell_metadata.txt.gz"

opts$classes <- c("WT", "Dnmt3a_KO", "Dnmt3b_KO", "Dnmt1_KO")

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  setnames("celltype.mapped","celltype") %>%
  .[pass_rnaQC==TRUE & class%in%opts$classes] %>%
  .[,c("cell", "class", "alias", "dataset", "celltype","nFeature_RNA","mit_percent_RNA", "rib_percent_RNA")] %>%
  setnames("alias","sample")
  # .[,celltype:=stringr::str_replace_all(celltype,opts$rename_celltypes)] %>%

table(sample_metadata$class)

##########
## Save ##
##########

fwrite(sample_metadata, io$outfile, sep="\t", na="NA", quote=F)
