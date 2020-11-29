library(ggpubr)
library(SingleCellExperiment)
library(scran)

#####################
## Define settings ##
#####################

# Load default settings
source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")

# I/O
io$outdir <- paste0(io$basedir,"/results/qc")
# io$seurat <- paste0(io$basedir,"/processed/all_batches/seurat.rds")
# io$metadata <- paste0(io$basedir,"/processed/all_batches/metadata.txt.gz")

# Options
opts$nFeature_RNA <- 2000
opts$nCount_RNA <- 4000
opts$percent.mt <- 10

###############
## Load data ##
###############

# Load Seurat object
srat <- readRDS(io$seurat)

# Load cell metadata
metadata <- fread(io$metadata)

#####################
## Plot QC metrics ##
#####################

to.plot <- srat@meta.data %>% as.data.table %>%
    melt(id.vars=c("batch","cell"), measure.vars=c("nCount_RNA","nFeature_RNA","percent.mt"))

tmp <- data.table(
    variable = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
    value = c(opts$nCount_RNA, opts$nFeature_RNA, opts$percent.mt)
)

p <- gghistogram(to.plot, x="value", fill="batch", bins=50) +
    geom_vline(aes(xintercept=value), linetype="dashed", data=tmp) +
    facet_wrap(~variable, scales="free") +
    theme(
        axis.text =  element_text(size=rel(0.8)),
        axis.title.x = element_blank(),
        legend.position = "right"
    )
    
pdf(sprintf("%s/qc_metrics.pdf",io$outdir), width=12, height=5, useDingbats = F)
print(p)
dev.off()

#############################
## Calculate doublet score ##
#############################

# opts$max_doublet_score <- 10000
# 
# sce <- as.SingleCellExperiment(srat)
# metadata[,doublet_score:=doubletCells(sce)]
# 
# # Plot
# to.plot <- metadata %>% 
#     melt(id.vars=c("batch","cell"), measure.vars=c("doublet_score"))
# 
# p <- gghistogram(to.plot, x="value", fill="batch", bins=50) +
#     geom_vline(xintercept=opts$max_doublet_score, linetype="dashed") +
#     theme(
#         axis.text =  element_text(size=rel(0.8)),
#         legend.position = "right"
#     )
# 
# table(metadata$doublet_score<opts$max_doublet_score)
# 
# pdf(sprintf("%s/doublet_score.pdf",io$outdir), width=9, height=5, useDingbats = F)
# print(p)
# dev.off()

############
## Filter ##
############

cells <- metadata %>%
    .[ nFeature_RNA>opts$nFeature_RNA & nCount_RNA>opts$nCount_RNA & percent.mt<opts$percent.mt,cell]
length(cells) / nrow(srat@meta.data)

# Subset cells that pass quality control
metadata[,pass_QC:=cell%in%cells]
srat <- srat[,cells]

##########
## Save ##
##########

saveRDS(srat, io$seurat)
fwrite(metadata, io$metadata, quote=F, na="NA", sep="\t")

