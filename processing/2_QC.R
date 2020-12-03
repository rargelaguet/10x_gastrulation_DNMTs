suppressPackageStartupMessages(library(ggpubr))

#####################
## Define settings ##
#####################

# Load default settings
source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")

# I/O
io$outdir <- paste0(io$basedir,"/results/qc/sixth_batch")
# io$seurat <- paste0(io$basedir,"/processed/all_batches/seurat.rds")
# io$metadata <- paste0(io$basedir,"/processed/all_batches/metadata.txt.gz")

# Options
opts$nFeature_RNA <- 1000
opts$nCount_RNA <- 2000
opts$percent.mt <- 15

###############
## Load data ##
###############

# io$metadata <- "/Users/ricard/data/10x_gastrulation_DNMTs/processed/all_batches/metadata.txt.gz"
io$metadata <- "/Users/ricard/data/10x_gastrulation_DNMTs/processed/sixth_batch/metadata.txt.gz"
metadata <- fread(io$metadata) %>% 
    # .[batch%in%opts$batches] %>%
    .[,pass_QC:=nFeature_RNA>opts$nFeature_RNA & nCount_RNA>opts$nCount_RNA & percent.mt<opts$percent.mt]

#####################
## Plot QC metrics ##
#####################

to.plot <- metadata %>%
    .[,log_nCount_RNA:=log2(nCount_RNA)] %>%
    melt(id.vars=c("batch","cell"), measure.vars=c("log_nCount_RNA","nFeature_RNA","percent.mt"))

## Box plot 

p <- ggboxplot(to.plot, x="batch", y="value") +
    facet_wrap(~variable, scales="free_y", nrow=3) +
    theme(
        axis.text.x = element_text(colour="black",size=rel(0.7), angle=90, hjust=1, vjust=0.5),  
        axis.title.x = element_blank()
    )

pdf(sprintf("%s/qc_metrics_boxplot.pdf",io$outdir), width=6, height=12, useDingbats = F)
print(p)
dev.off()

## histogram 

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



########################################################
## Plot fraction of cells that pass QC for each batch ##
########################################################

to.plot <- metadata %>%
    .[,mean(pass_QC),by="batch"]

p <- ggbarplot(to.plot, x="batch", y="V1", fill="gray70") +
    labs(x="", y="Fraction of cells that pass QC") +
    # facet_wrap(~stage)
    theme(
        axis.text.x = element_text(colour="black",size=rel(0.75), angle=90, hjust=1, vjust=0.5)
    )

pdf(sprintf("%s/qc_metrics_barplot.pdf",io$outdir), width=9, height=7, useDingbats = F)
print(p)
dev.off()

##########
## Save ##
##########

fwrite(metadata, io$metadata, quote=F, na="NA", sep="\t")

