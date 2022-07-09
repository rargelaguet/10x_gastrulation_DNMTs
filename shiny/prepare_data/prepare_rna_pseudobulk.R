source(here::here("settings.R"))
source(here::here("utils.R"))
source(here::here("pseudobulk/utils.R"))

# I/O
io$metadata <- file.path(io$basedir,"results/mapping/sample_metadata_after_mapping.txt.gz")
io$sce <- file.path(io$basedir,"processed/SingleCellExperiment.rds")
io$outfile <- "/Users/argelagr/data/shiny_dnmt_tet/SingleCellExperiment_shiny.rds"


# Options
opts$group_by <- "class_sample_celltype_dataset"

opts$rename_celltypes <- c(
  "Erythroid3" = "Erythroid",
  "Erythroid2" = "Erythroid",
  "Erythroid1" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors",
  "Anterior_Primitive_Streak" = "Primitive_Streak"
)

opts$classes <- c(
  "WT",
  "Dnmt3a_KO",
  "Dnmt3b_KO",
  "Dnmt1_KO"
  # "Dnmt3ab_KO"
)


########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  setnames("celltype.mapped","celltype") %>%
  .[,dataset:=ifelse(grepl("Grosswendt",sample),"CRISPR","KO")] %>%
  .[,celltype:=stringr::str_replace_all(celltype,opts$rename_celltypes)] %>%
  .[,class_celltype:=sprintf("%s-%s",class,celltype)] %>%
  .[,class_celltype_dataset:=sprintf("%s-%s-%s",class,celltype,dataset)] %>%
  .[,class_sample_celltype_dataset:=sprintf("%s-%s-%s-%s",class,sample,celltype,dataset)] %>%
  # .[,class_celltype_dataset:=sprintf("%s-%s-%s",class,celltype,dataset)] %>%
  .[pass_rnaQC==TRUE & !is.na(eval(as.name(opts$group_by)))]

######################################################
## Calculate pseudovulk stats and do some filtering ##
######################################################

pseudobulk_stats.dt <- sample_metadata[,.N,by=c("class","sample","celltype","dataset","class_sample_celltype_dataset")]

# Filter each instance by minimum number of cells
pseudobulk_stats.dt <- pseudobulk_stats.dt[N>=30]

# Select celltypes that are measured in at least 5 WT samples
celltypes.to.use <- pseudobulk_stats.dt[class=="WT",.N,by=c("celltype")] %>% .[N>=5,celltype] 
pseudobulk_stats.dt <- pseudobulk_stats.dt[celltype%in%celltypes.to.use]

# For each class and celltype combination, require at least 3 samples
tmp <- pseudobulk_stats.dt[,.N,by=c("celltype","class")] %>% .[N>=3] %>% .[,N:=NULL]
pseudobulk_stats.dt <- pseudobulk_stats.dt %>% merge(tmp,by=c("celltype","class"))

# Print stats
# print(pseudobulk_stats.dt)

# Update metadata
sample_metadata <- sample_metadata[class_sample_celltype_dataset%in%pseudobulk_stats.dt$class_sample_celltype_dataset]

##############################
## Load RNA expression data ##
##############################

sce <- load_SingleCellExperiment(io$sce, cells=sample_metadata$cell)
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

################
## Pseudobulk ##
################

sce_pseudobulk <- pseudobulk_sce_fn(
  x = sce,
  assay = "counts",
  by = opts$group_by,
  fun = "sum",
  scale = FALSE # Should pseudo-bulks be scaled with the effective library size & multiplied by 1M?
)

assayNames(sce_pseudobulk) <- "counts"

###############
## Parse SCE ##
###############

# Add metadata
sce_pseudobulk$class <- stringr::str_split(colnames(sce_pseudobulk), pattern = "-") %>% map_chr(1)
sce_pseudobulk$sample <- stringr::str_split(colnames(sce_pseudobulk), pattern = "-") %>% map_chr(2)
sce_pseudobulk$celltype <- stringr::str_split(colnames(sce_pseudobulk), pattern = "-") %>% map_chr(3)
sce_pseudobulk$dataset <- stringr::str_split(colnames(sce_pseudobulk), pattern = "-") %>% map_chr(4)

# Filter genes
genes <- fread("/Users/argelagr/data/shiny_dnmt_tet/genes.txt", header = F)[[1]]
sce_pseudobulk <- sce_pseudobulk[genes,]

# Sanity checks
stopifnot(unique(sce_pseudobulk$class)%in%opts$classes)
stopifnot(unique(sce_pseudobulk$sample)%in%opts$samples)
stopifnot(unique(sce_pseudobulk$celltype)%in%c(opts$celltypes,c("Erythroid", "Blood_progenitors")))
stopifnot(unique(sce_pseudobulk$dataset)%in%c("KO","CRISPR"))

###################
## Normalisation ##
###################

logcounts(sce_pseudobulk) <- log2(1e6*(sweep(counts(sce_pseudobulk),2,colSums(counts(sce_pseudobulk)),"/"))+1)

# Remove counts assay
assays(sce_pseudobulk)["counts"] <- NULL

# Save
saveRDS(sce_pseudobulk, io$outfile)

#####################
## Save statistics ##
#####################

# stats.dt <- data.table(table(sce[[opts$group_by]])) %>% setnames(c("sample","N"))
# fwrite(stats.dt, file.path(io$outdir,sprintf("pseudobulk_stats_%s.txt.gz",opts$group_by)), sep="\t")
