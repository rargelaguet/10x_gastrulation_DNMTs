################################################################################
## Sex determination using RNA expression counts of genes on the Y chromosome ##
################################################################################

suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(argparse))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--samples',   type="character",   nargs='+',  help='samples')
p$add_argument('--sce',       type="character",               help='SingleCellExperiment file')
p$add_argument('--metadata',  type="character",               help='metadata file')
p$add_argument('--threshold.ratioY',  type="double", default=1e-3,              help='XXX')
p$add_argument('--outdir',          type="character",               help='Output directory')
p$add_argument('--test',            action = "store_true",          help='Testing mode')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
} else {
  source("/homes/ricard/10x_gastrulation_DNMTs/settings.R")
}

## START TEST ##
# args$samples <- opts$batches
# args$sce <- io$sce
# args$metadata <- io$metadata
# args$outdir <- paste0(io$basedir,"/results/sex")
# args$test <- TRUE
# args$threshold.ratioY <- 1e-3
## END TEST ##

####################
## Define options ##
####################


###############
## Load data ##
###############

if (isTRUE(args$test)) args$samples <- head(args$samples,n=2)

# Load sample metadata
sample_metadata <- fread(args$metadata) %>% .[pass_QC==T & batch%in%args$samples]

# Load SingleCellExperiment
# sce <- readRDS(io$sce)[,sample_metadata$cell]
sce <- load_SingleCellExperiment(args$sce, cells = sample_metadata$cell)
dim(sce)

# Load gene metadata
gene_metadata <- fread(io$gene_metadata) %>% 
  .[symbol%in%rownames(sce)]

################
## Parse data ##
################

# Group genes by chromosome
genes.chrY <- gene_metadata[chr=="chrY",symbol]
genes.chrX <- gene_metadata[chr=="chrX",symbol]
genes.chr10 <- gene_metadata[chr=="chr10",symbol]

# Manual filtering
# For some reason Erdr1 is predicted as Ychr, but the last version of ENSEMBL is in the Xchr
genes.chrY <- genes.chrY[!genes.chrY=="ENSMUSG00000096768"]


# Create data.table from SingleCellExperiment object
# dt <- counts(sce[c(genes.chrY,genes.chrX,genes.chr10)]) %>% as.matrix %>% t %>%
#   as.data.table(keep.rownames="cell") %>% 
#   melt(id.vars="cell", value.name="counts", variable.name="symbol") %>%
#   merge(sample_metadata[,c("cell","batch")], by="cell")#  %>%
#   # merge(gene_metadata[,c("chr","ens_id","symbol")], by="symbol")
dt <- args$samples %>% map(function(i) {
  sce.filt <- sce[,sce$batch==i] %>% .[c(genes.chrY,genes.chrX,genes.chr10),]
  dt <- data.table(
    symbol = rownames(sce.filt),
    counts = counts(sce.filt) %>% Matrix::rowSums()
  ) %>% .[,batch:=i]
}) %>% rbindlist %>% merge(gene_metadata[,c("chr","ens_id","symbol")], by="symbol")

to.plot <- dt[chr=="chrY"] %>% 
  # .[,.(counts=sum(counts)),by=c("batch","chr","ens_id","symbol")] %>%
  .[,.(counts=sum(counts)),by=c("batch","chr","symbol")] %>%
  .[,mean:=mean(counts),by="symbol"] %>% .[mean>0] %>% .[,mean:=NULL] 

##########
## Plot ##
##########

p <- ggbarplot(to.plot, x="symbol", y="counts", facet="batch", fill="gray70") +
# p <- ggbarplot(to.plot, x="ens_id", y="counts", facet="batch", fill="gray70") +
  labs(x="", y="Read counts") +
  theme(
  axis.text.x = element_text(colour="black",size=rel(0.7), angle=45, hjust=1, vjust=1),
  axis.text.y = element_text(colour="black",size=rel(0.8))
  )

pdf(sprintf("%s/pdf/sex_ychr_expr_per_gene.pdf",args$outdir))
print(p)
dev.off()

###################################################
## Barplots of chrXY/chr1 count ratio per embryo ##
###################################################

# Agregate counts over all genes
to.plot <- dt %>% 
  .[,.(counts=sum(counts)),by=c("batch","chr")] %>%
  dcast(batch~chr, value.var="counts") %>%
  .[,ratioY:=chrY/chr10] %>% .[,ratioX:=chrX/chr10] %>%
  .[,sex:=c("female","male")[as.numeric(ratioY>=args$threshold.ratioY)+1]]

p <- ggbarplot(to.plot, x="batch", y="ratioY", fill="sex", sort.val = "asc", palette="Dark2") +
  labs(x="", y="chrY/chr1 counts ratio") +
  theme(
    legend.position = "right",
    # axis.text.x = element_text(colour="black",size=rel(0.8), angle=45, hjust=1, vjust=1),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

pdf(sprintf("%s/pdf/sex_ychr_expr_aggregated.pdf",args$outdir))
print(p)
dev.off()


#############################
## Plot expression of Xist ##
#############################

to.plot <- dt[symbol=="Xist"] %>% 
  merge(to.plot[,c("sex","batch")], by="batch")

p <- ggbarplot(to.plot, x="batch", y="counts", fill="sex") +
  labs(x="", y="Xist expression") +
  theme(
    axis.text.x = element_text(colour="black",size=rel(0.8), angle=45, hjust=1, vjust=1),
  )

pdf(sprintf("%s/pdf/Xist_expr.pdf",args$outdir))
print(p)
dev.off()

##################################################################
## Scatterplot of chrY/chr1 vs chrX/chr1 count ratio per embryo ##
##################################################################

# ggscatter(to.plot, x="ratioX", y="ratioY") +
#   labs(x="chrX/chr1 ratio", y="chrY/chr1 ratio") +
#   theme(
#     # axis.text.x = element_text(colour="black",size=rel(0.8), angle=45, hjust=1, vjust=1),
#     # axis.title.x = element_blank(),
#     # axis.ticks.x = element_blank()
#   )

##########
## Save ##
##########

fwrite(to.plot, paste0(args$outdir,"/sex_assignment.txt.gz"))





