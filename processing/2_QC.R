suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(argparse))

#####################
## Define arguments ##
#####################

p <- ArgumentParser(description='')
p$add_argument('--outputdir',       type="character",                    help='Output directory')
p$add_argument('--nFeature_RNA',       type="integer",                    help='Minimum number of expressed genes')
p$add_argument('--nCount_RNA',       type="integer",                    help='Minimum number of reads')
p$add_argument('--percent.mt',       type="integer",                    help='Maximum percentage of mitochondrial reads')
p$add_argument('--samples',         type="character",       nargs="+",   help='Samples')
args <- p$parse_args(commandArgs(TRUE))


#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
    source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
    source("/homes/ricard/10x_gastrulation_DNMTs/settings.R")
}

## START TEST ##
# args <- list()
# args$outputdir <- paste0(io$basedir,"/results/qc")
# args$samples <- head(opts$batches,n=2)
# args$nFeature_RNA <- 1000
# args$nCount_RNA <- 2000
# args$percent.mt <- 15
## END TEST ##

# Sanity checks
stopifnot(args$samples%in%opts$batches)

###############
## Load data ##
###############

# io$metadata <- "/Users/ricard/data/10x_gastrulation_DNMTs/processed/all_batches/metadata.txt.gz"
# io$metadata <- "/Users/ricard/data/10x_gastrulation_DNMTs/processed/sixth_batch/metadata.txt.gz"
metadata <- fread(io$metadata) %>% 
    .[batch%in%args$samples] %>%
    .[,pass_QC:=nFeature_RNA>args$nFeature_RNA & nCount_RNA>args$nCount_RNA & percent.mt<args$percent.mt]

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
        axis.text.x = element_text(colour="black",size=rel(0.65), angle=20, hjust=1, vjust=1),  
        axis.title.x = element_blank()
    )

# pdf(sprintf("%s/qc_metrics_boxplot.pdf",args$outputdir), width=6, height=12, useDingbats = F)
pdf(sprintf("%s/qc_metrics_boxplot.pdf",args$outputdir), useDingbats = F)
print(p)
dev.off()

## histogram 

tmp <- data.table(
    variable = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
    value = c(args$nCount_RNA, args$nFeature_RNA, args$percent.mt)
)

p <- gghistogram(to.plot, x="value", fill="batch", bins=50) +
    geom_vline(aes(xintercept=value), linetype="dashed", data=tmp) +
    facet_wrap(~variable, scales="free") +
    theme(
        axis.text =  element_text(size=rel(0.8)),
        axis.title.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(size=rel(0.5))
    )
    
# pdf(sprintf("%s/qc_metrics.pdf",args$outputdir), width=12, height=5, useDingbats = F)
pdf(sprintf("%s/qc_metrics_histogram.pdf",args$outputdir), useDingbats = F)
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
        axis.text.x = element_text(colour="black",size=rel(0.65), angle=20, hjust=1, vjust=1),  
    )

# pdf(sprintf("%s/qc_metrics_barplot.pdf",args$outputdir), width=9, height=7, useDingbats = F)
pdf(sprintf("%s/qc_metrics_barplot.pdf",args$outputdir), useDingbats = F)
print(p)
dev.off()

##########
## Save ##
##########

fwrite(metadata, paste0(args$outputdir,"/sample_metadata_after_qc.txt.gz"), quote=F, na="NA", sep="\t")

