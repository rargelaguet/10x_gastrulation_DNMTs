here::i_am("mapping/run/parse_sample_metadata_after_mapping.R")

source(here::here("settings.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',    type="character",  help='Metadata file to use as input')
# p$add_argument('--mapping_seurat',    type="character", nargs="+", help='Results of the Seurat mapping')
p$add_argument('--mapping_mnn',    type="character",  nargs="+", help='Results of the MNN mapping')
p$add_argument('--outfile',          type="character",               help='Output file')
args <- p$parse_args(commandArgs(TRUE))

###################
## Load settings ##
###################

## START TEST ##
# args$metadata <- file.path(io$basedir,"results_new/doublet_detection/sample_metadata_after_doublets.txt.gz")
# args$mapping_mnn <- file.path(io$basedir,"results_new/mapping/mapping_mnn_all_samples.txt.gz")
# args$outfile <- file.path(io$basedir,"results_new/mapping/sample_metadata_after_mapping.txt.gz")
## END TEST ##

###################
## Load metadata ##
###################

sample_metadata <- fread(args$metadata)

# Temporary
stopifnot(sample_metadata$sample%in%names(opts$sample2class))
sample_metadata[,class:=stringr::str_replace_all(sample,opts$sample2class)]
stopifnot(!is.na(sample_metadata$class))
print(table(sample_metadata$class))

##########################
## Load mapping results ##
##########################

# MNN
mapping_mnn.dt <- args$mapping_mnn %>% map(~ fread(.)) %>% rbindlist
stopifnot(mapping_mnn.dt$cell%in%sample_metadata$cell)

# Seurat
# mapping_seurat.dt <- args$mapping_seurat %>% map(~ fread(.)) %>% rbindlist
# stopifnot(mapping_seurat.dt$cell%in%sample_metadata$cell)

###########
## Merge ##
###########

# mapping.dt <- merge(mapping_mnn.dt, mapping_seurat.dt, by="cell", suffixes=c("_mnn","_seurat"))
# to.save <- sample_metadata %>% merge(mapping.dt, by="cell", all.x=TRUE)

to.save <- sample_metadata %>% merge(mapping_mnn.dt, by="cell", all.x=TRUE)

#######################
## Rename cell types ##
#######################

# opts$rename.celltypes <- c(
#   "Erythroid1" = "early_Erythroid",
#   "Erythroid2" = "early_Erythroid",
#   "Erythroid3" = "late_Erythroid",
#   "Blood_progenitors_1" = "Blood_progenitors",
#   "Blood_progenitors_2" = "Blood_progenitors",
#   "Rostral_neurectoderm" = "Neurectoderm",
#   "Caudal_neurectoderm" = "Neurectoderm",
#   "Anterior_Primitive_Streak" = "Primitive_Streak",
#   "Mixed_mesoderm" = "Nascent_mesoderm",
#   "Allantois" = "ExE_mesoderm"
# )

# to.save %>% .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,opts$rename.celltypes)]

to.save %>%
  .[,celltype_class:=as.character(NA)] %>%
  .[!is.na(celltype.mapped),celltype_class:=sprintf("%s_%s",celltype.mapped,class)]

# stopifnot(!is.na(to.save$celltype_class))
# print(table(to.save$celltype_class))

#################
## Save output ##
#################

fwrite(to.save, args$outfile, sep="\t", na="NA", quote=F)

######################
## Compare mappings ##
######################

# mapping_mnn.dt <- readRDS(sprintf("%s/mapping_mnn_%s.rds",io$mapping.dir,paste(opts$samples,collapse="-")))$mapping %>% .[,c("cell","celltype.mapped","celltype.score","closest.cell")] %>% as.data.table
# mapping_seurat.dt <- fread(sprintf("%s/mapping_seurat_%s.txt.gz",io$mapping.dir,paste(opts$samples,collapse="-"))) %>% .[,c("predicted.id")] %>% as.data.table
# 
# foo <- merge(
#   mapping_mnn.dt[,c("cell","celltype.mapped")] %>% setnames("celltype.mapped","celltype_mnn"),
#   mapping_seurat.dt[,c("cell","predicted.id")] %>% setnames("predicted.id","celltype_seurat"),
#   by = c("cell")
# )
