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
  source("/Users/ricard/10x_gastrulation_DNMTs/utils.R")
} else {
  source("/homes/ricard/10x_gastrulation_DNMTs/settings.R")
  source("/homes/ricard/10x_gastrulation_DNMTs/utils.R")
}

## START TEST ##
args$samples <- opts$samples
args$sce <- io$sce
args$metadata <- io$metadata
args$outdir <- paste0(io$basedir,"/results/sex")
args$test <- FALSE
args$threshold.ratioY <- 1e-3
## END TEST ##

####################
## Define options ##
####################


###############
## Load data ##
###############

if (isTRUE(args$test)) args$samples <- head(args$samples,n=2)

# Load sample metadata
sample_metadata <- fread(args$metadata) %>% .[pass_QC==T & sample%in%args$samples]

# Load SingleCellExperiment
# sce <- readRDS(io$sce)[,sample_metadata$cell]
sce <- load_SingleCellExperiment(args$sce, cells = sample_metadata$cell)
dim(sce)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

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
#   merge(sample_metadata[,c("cell","sample")], by="cell")#  %>%
#   # merge(gene_metadata[,c("chr","ens_id","symbol")], by="symbol")
dt <- args$samples %>% map(function(i) {
  sce.filt <- sce[,sce$sample==i] %>% .[c(genes.chrY,genes.chrX,genes.chr10),]
  dt <- data.table(
    symbol = rownames(sce.filt),
    counts = counts(sce.filt) %>% Matrix::rowSums()
  ) %>% .[,sample:=i]
}) %>% rbindlist %>% merge(gene_metadata[,c("chr","ens_id","symbol")], by="symbol")

to.plot <- dt[chr=="chrY"] %>% 
  # .[,.(counts=sum(counts)),by=c("sample","chr","ens_id","symbol")] %>%
  .[,.(counts=sum(counts)),by=c("sample","chr","symbol")] %>%
  .[,mean:=mean(counts),by="symbol"] %>% .[mean>0] %>% .[,mean:=NULL] 

##########
## Plot ##
##########

p <- ggbarplot(to.plot, x="symbol", y="counts", facet="sample", fill="gray70") +
# p <- ggbarplot(to.plot, x="ens_id", y="counts", facet="sample", fill="gray70") +
  labs(x="", y="Read counts") +
  guides(x = guide_axis(angle = 90)) +
  theme(
  axis.text.x = element_text(colour="black",size=rel(0.5)),
  axis.ticks.x = element_line(size=rel(0.5)),
  axis.text.y = element_text(colour="black",size=rel(0.8))
  )

pdf(sprintf("%s/pdf/sex_ychr_expr_per_gene.pdf",args$outdir), width=16, height=10)
print(p)
dev.off()

###################################################
## Barplots of chrXY/chr1 count ratio per embryo ##
###################################################

# Agregate counts over all genes
sex_assignment.dt <- dt %>% 
  .[,.(counts=sum(counts)),by=c("sample","chr")] %>%
  dcast(sample~chr, value.var="counts") %>%
  .[,ratioY:=chrY/chr10] %>% .[,ratioX:=chrX/chr10] %>%
  .[,sex:=c("female","male")[as.numeric(ratioY>=args$threshold.ratioY)+1]]

p <- ggbarplot(sex_assignment.dt, x="sample", y="ratioY", fill="sex", sort.val = "asc", palette="Dark2") +
  labs(x="", y="chrY/chr1 counts ratio") +
  guides(x = guide_axis(angle = 90)) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(colour="black",size=rel(0.6))
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank()
  )

pdf(sprintf("%s/pdf/sex_ychr_expr_aggregated.pdf",args$outdir))
print(p)
dev.off()


#############################
## Plot expression of Xist ##
#############################

to.plot <- dt[symbol=="Xist"] %>% 
  # .[,expr:=log(counts+1)] %>%
  merge(sex_assignment.dt[,c("sex","sample")], by="sample")

p <- ggbarplot(to.plot, x="sample", y="counts", fill="sex") +
  labs(x="", y="Xist expression") +
  guides(x = guide_axis(angle = 90)) +
  theme(
    axis.text.x = element_text(colour="black",size=rel(0.6))
    # axis.text.x = element_text(colour="black",size=rel(0.8), angle=45, hjust=1, vjust=1),
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

############################
## Update sample metadata ##
############################

sample_metadata <- fread(io$metadata)
sample_metadata <- sample_metadata %>% merge(sex_assignment.dt[,c("sample","sex")], by="sample", all.x = TRUE)
fwrite(sample_metadata, io$metadata, sep="\t", quote=F, na="NA")

##########
## Save ##
##########

fwrite(sex_assignment.dt[,c("sample","sex")], paste0(args$outdir,"/sex_assignment.txt.gz"))





