suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(SingleCellExperiment))

################
## Define I/O ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
  source("/Users/ricard/10x_gastrulation_DNMTs/mapping/run/mapping_functions.R")
} else {
  source("/homes/ricard/10x_gastrulation_DNMTs/settings.R")
  source("/homes/ricard/10x_gastrulation_DNMTs/mapping/run/mapping_functions.R")
}
io$path2atlas <- io$atlas.basedir
io$path2query <- io$basedir
io$outdir <- paste0(io$basedir,"/results/second_batch/mapping")

####################
## Define options ##
####################

opts$atlas.stages <- c(
  "E6.5",
  "E6.75",
  "E7.0",
  "E7.25",
  "E7.5",
  "E7.75",
  "E8.0",
  "E8.25",
  "E8.5",
  "mixed_gastrulation"
)

opts$query.batches <- c(
  # "E75_TET_TKO_L002", 
  # "E75_WT_Host_L001", 
  # "E85_Rep1_TET_TKO_L004", 
  # "E85_Rep2_TET_TKO_L006", 
  "E85_Rep1_WT_Host_L003"
  # "E85_Rep2_WT_Host_L005"
  # "E125_DNMT3A_HET_A_L001", 
  # "E125_DNMT3A_HET_A_L003",
  # "E125_DNMT3A_KO_B_L002", 
  # "E125_DNMT3A_KO_E_L004"
)

################
## Load atlas ##
################

sce_atlas  <- readRDS(io$atlas.sce)
meta_atlas <- fread(io$atlas.metadata) %>%
  .[stripped==F & doublet==F]

# Filter
meta_atlas <- meta_atlas[meta_atlas$stage%in%opts$atlas.stages,]
sce_atlas <- sce_atlas[,meta_atlas$cell] 

################
## Load query ##
################

# Load Seurat object
sce_query <- readRDS(io$sce)
meta_query <- fread(io$metadata) %>%
  .[pass_QC==T & batch%in%opts$query.batches & cell%in%colnames(sce_query)]


# Filter
sce_query <- sce_query[,meta_query$cell]

#############
## Prepare ## 
#############

genes.intersect <- intersect(rownames(sce_query), rownames(sce_atlas))
sce_query  <- sce_query[genes.intersect,]
sce_atlas <- sce_atlas[genes.intersect,]

#########
## Map ##
#########

marker_genes.dt <- fread(io$atlas.marker_genes)
marker_genes.dt <- marker_genes.dt[,head(.SD,n=50),by="celltype"]
marker_genes <- marker_genes.dt$ens_id
marker_genes <- marker_genes[marker_genes%in%genes.intersect]

print(marker_genes.dt[,.N,by="celltype"])

mapping  <- mapWrap(
  atlas_sce = sce_atlas, 
  atlas_meta = meta_atlas,
  map_sce = sce_query, 
  map_meta = meta_query, 
  genes = marker_genes,
  order = NULL,
  npcs = 50,
  k = 25,
  return.list = FALSE
)

##########
## Save ##
##########

# save mapping results as an .rds file
saveRDS(mapping, paste0(io$outdir,"/mapping10x_mnn.rds"))
