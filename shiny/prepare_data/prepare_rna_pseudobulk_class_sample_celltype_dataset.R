source(here::here("settings.R"))
source(here::here("utils.R"))
source(here::here("pseudobulk/utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$metadata <- file.path(io$basedir,"results/mapping/sample_metadata_after_mapping.txt.gz")
io$sce <- file.path(io$basedir,"processed/SingleCellExperiment.rds")
io$genes <- file.path(io$basedir,"shiny/rna_expression/genes.txt")
io$outdir <- file.path(io$basedir,"shiny/rna_expression/pseudobulk/class_sample_celltype_dataset"); dir.create(io$outdir, showWarnings = F)

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

cell_metadata.dt <- fread(io$metadata) %>%
  setnames("celltype.mapped","celltype") %>%
  .[,dataset:=ifelse(grepl("Grosswendt",sample),"CRISPR","KO")] %>%
  .[,celltype:=stringr::str_replace_all(celltype,opts$rename_celltypes)] %>%
  .[,class_sample_celltype_dataset:=sprintf("%s-%s-%s-%s",class,sample,celltype,dataset)] %>%
  .[pass_rnaQC==TRUE & !is.na(eval(as.name(opts$group_by)))]

######################################################
## Calculate pseudovulk stats and do some filtering ##
######################################################

pseudobulk_stats.dt <- cell_metadata.dt[,.N,by=c("class","sample","celltype","dataset","class_sample_celltype_dataset")]

# Filter each instance by minimum number of cells
pseudobulk_stats.dt <- pseudobulk_stats.dt[N>=25]

# Select celltypes that are measured in at least 5 WT samples
celltypes.to.use <- pseudobulk_stats.dt[class=="WT",.N,by=c("celltype")] %>% .[N>=5,celltype] 
pseudobulk_stats.dt <- pseudobulk_stats.dt[celltype%in%celltypes.to.use]

# For each class and celltype combination, require at least 3 samples
tmp <- pseudobulk_stats.dt[,.N,by=c("celltype","class")] %>% .[N>=3] %>% .[,N:=NULL]
pseudobulk_stats.dt <- pseudobulk_stats.dt %>% merge(tmp,by=c("celltype","class"))

# Print stats
# print(pseudobulk_stats.dt)

# Update metadata
cell_metadata.dt <- cell_metadata.dt[class_sample_celltype_dataset%in%pseudobulk_stats.dt$class_sample_celltype_dataset]

##############################
## Load RNA expression data ##
##############################

sce <- load_SingleCellExperiment(io$sce, cells=cell_metadata.dt$cell)
colData(sce) <- cell_metadata.dt %>% tibble::column_to_rownames("cell") %>% DataFrame

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

##################
## Subset genes ##
##################

genes <- fread(io$genes, header = F)[[1]]
sce_pseudobulk <- sce_pseudobulk[genes,]

###############
## Parse SCE ##
###############

# # Add metadata
# sce_pseudobulk$class <- stringr::str_split(colnames(sce_pseudobulk), pattern = "-") %>% map_chr(1)
# sce_pseudobulk$sample <- stringr::str_split(colnames(sce_pseudobulk), pattern = "-") %>% map_chr(2)
# sce_pseudobulk$celltype <- stringr::str_split(colnames(sce_pseudobulk), pattern = "-") %>% map_chr(3)
# sce_pseudobulk$dataset <- stringr::str_split(colnames(sce_pseudobulk), pattern = "-") %>% map_chr(4)
# 
# # Sanity checks
# stopifnot(unique(sce_pseudobulk$class)%in%opts$classes)
# stopifnot(unique(sce_pseudobulk$sample)%in%opts$samples)
# stopifnot(unique(sce_pseudobulk$celltype)%in%c(opts$celltypes,c("Erythroid", "Blood_progenitors")))
# stopifnot(unique(sce_pseudobulk$dataset)%in%c("KO","CRISPR"))

###################
## Normalisation ##
###################

logcounts(sce_pseudobulk) <- log2(1e6*(sweep(counts(sce_pseudobulk),2,colSums(counts(sce_pseudobulk)),"/"))+1)

# Remove counts assay
# assays(sce_pseudobulk)["counts"] <- NULL

##########
## Save ##
##########

# Save sample metadata
to.save <- pseudobulk_stats.dt %>% copy %>% setnames("class_sample_celltype_dataset","id")
fwrite(to.save, file.path(io$outdir,"sample_metadata.txt.gz"), na="NA", quote=F, sep="\t")
       
# Save expression matrix
saveRDS(round(logcounts(sce_pseudobulk),2), file.path(io$outdir,"rna_expr.rds"))
