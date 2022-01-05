here::i_am("processing/1_create_seurat_rna.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(Seurat))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--inputdir',        type="character",                    help='Input directory')
p$add_argument('--outdir',       type="character",                    help='Output directory')
p$add_argument('--samples',         type="character",       nargs="+",   help='Samples')
p$add_argument('--test',            action="store_true",                 help='Testing mode')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
# args <- list()
# args$inputdir <- paste0(io$basedir,"/original")
# args$outdir <- paste0(io$basedir,"/processed_merged")
# args$samples <- c("Dnmt3a_E8.5_embryo3_Grosswendt2020","WT_E8.5_embryo9_Grosswendt2020") # opts$samples[1:2]
# args$test <- FALSE
## END TEST ##


# These settings are important for consistency with ArchR, which provides little flexibility to edit cell names
opts$trim.barcode <- FALSE
opts$sample_cell_separator <- "_"
opts$min.counts <- 500

##############################
## Load and merge data sets ##
##############################

stopifnot(args$samples%in%opts$samples)
if (args$test) args$samples <- head(args$samples,n=2)

count_mtx <- list()
cell.info <- list()

for (i in args$samples) {
  print(i)
    
  # Load cell metadata
  barcode.loc <- file.path(args$inputdir,sprintf("%s/barcodes.tsv.gz",i))
  cell.info[[i]] <- fread(barcode.loc, header=F) %>%
    setnames("barcode") %>%
    .[,barcode:=ifelse(rep(opts$trim.barcode,.N),gsub("-1","",barcode),barcode)] %>%
    .[,c("sample","cell"):=list(i,sprintf("%s%s%s",i,opts$sample_cell_separator,barcode))]
  dim(cell.info[[i]])
  stopifnot(!duplicated(cell.info[[i]]$barcode))
  
  # Load matrix  
  count_mtx[[i]] <- Read10X(file.path(args$inputdir,i))
  
  # Do initial filtering of low-quality cells
  count_mtx[[i]] <- count_mtx[[i]][,colSums(count_mtx[[i]])>=opts$min.counts]
  cell.info[[i]] <- cell.info[[i]][barcode%in%colnames(count_mtx[[i]])] %>% setkey(barcode) %>% .[colnames(count_mtx[[i]])]
  stopifnot(complete.cases(cell.info[[i]]))
  stopifnot(colnames(count_mtx[[i]])==cell.info[[i]]$barcode)
  
  # remove human genes for samples that are mapped to the mixed transcriptome
  if (any(grep("^mm10_",rownames(count_mtx[[i]])))) {
    count_mtx[[i]] <- count_mtx[[i]][grep("mm10",rownames(count_mtx[[i]])),]
    rownames(count_mtx[[i]]) <- rownames(count_mtx[[i]]) %>% gsub("mm10___","",.)
  }
  
}

print(lapply(count_mtx,dim))

#######################
## Keep common genes ##
#######################

genes <- Reduce("intersect",lapply(count_mtx,rownames))
for (i in 1:length(count_mtx)) {
  count_mtx[[i]] <- count_mtx[[i]][genes,]
}

stopifnot(length(unique(lapply(count_mtx,nrow)))==1)
stopifnot(length(unique(lapply(count_mtx,rownames)))==1)

#################
## Concatenate ##
#################

# Concatenate cell metadata
cell.info <- rbindlist(cell.info)
rownames(cell.info) <- cell.info$cell

# Concatenate matrices
count_mtx <- do.call("cbind",count_mtx)
colnames(count_mtx) <- cell.info$cell

##################
## Filter genes ##
##################

# Remove duplicated genes
count_mtx <- count_mtx[!duplicated(rownames(count_mtx)),]

# Sanity checks
stopifnot(sum(duplicated(rownames(count_mtx)))==0)
stopifnot(sum(duplicated(colnames(count_mtx)))==0)

##########################
## Create Seurat object ##
##########################

cell.info.to.seurat <- cell.info[cell%in%colnames(count_mtx)] %>% setkey(cell) %>% .[colnames(count_mtx)] %>% as.data.frame
rownames(cell.info.to.seurat) <- cell.info.to.seurat$cell
stopifnot(rownames(cell.info.to.seurat)==colnames(count_mtx))
stopifnot(sum(is.na(rownames(cell.info.to.seurat$cell)))==0)

seurat <- CreateSeuratObject(count_mtx, meta.data = cell.info.to.seurat)

# Add mit percenatge
seurat[["mit_percent_RNA"]] <- PercentageFeatureSet(seurat, pattern = "^mt-") %>% round(2)

# Add rib RNA content
ribo.genes <- grep(pattern = "^Rp[l|s]", x = rownames(seurat), value = TRUE)
seurat[["rib_percent_RNA"]] <- PercentageFeatureSet(seurat, features = ribo.genes) %>% round(2)

#####################
## Create metadata ##
#####################

metadata <- seurat@meta.data %>% as.data.table %>% .[,orig.ident:=NULL] %>%
  .[,c("cell","barcode","sample","nFeature_RNA","nCount_RNA","mit_percent_RNA","rib_percent_RNA")]

# Add stage information
metadata[,stage:=as.character(NA)] %>%
  # .[grepl("E8",sample),stage:="E8.5"] %>%
  .[,stage:="E8.5"]# %>%
  # .[grepl("E12",sample),stage:="E12.5"]
print(table(metadata$stage))
stopifnot(!is.na(metadata$stage))

# Add class information
stopifnot(metadata$sample%in%names(opts$sample2class))
metadata[,class:=stringr::str_replace_all(sample,opts$sample2class)]
print(table(metadata$class))
stopifnot(!is.na(metadata$class))

# Add alias information
stopifnot(metadata$sample%in%names(opts$sample2alias))
metadata[,alias:=stringr::str_replace_all(sample,opts$sample2alias)]
print(table(metadata$alias))
stopifnot(!is.na(metadata$alias))

# Parse metadata Grosswendt2020
metadata[,dataset:=ifelse(grepl("Grosswendt",sample),"Grosswendt","This_data_set")]

##########
## Save ##
##########

fwrite(metadata, file.path(args$outdir,"metadata.txt.gz"), quote=F, na="NA", sep="\t")
# fwrite(cell.info, paste0(args$outdir,"/cell_info.txt.gz"), quote=F, na="NA", sep="\t")
# fwrite(gene.info, paste0(args$outdir,"/gene_info.txt.gz"), quote=F, na="NA", sep="\t")
saveRDS(seurat, file.path(args$outdir,"seurat.rds"))

