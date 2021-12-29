here::i_am("sex/sex_determination.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--samples',   type="character",   nargs='+',  help='samples')
p$add_argument('--sce',       type="character",               help='SingleCellExperiment file')
p$add_argument('--metadata',  type="character",               help='metadata file')
p$add_argument('--threshold.ratioY',  type="double", default=1e-3,              help='XXX')
p$add_argument('--outdir',          type="character",               help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
args$samples <- opts$samples[1:2]
args$sce <- io$sce
args$metadata <- io$metadata
args$outdir <- paste0(io$basedir,"/results_new/sex")
args$threshold.ratioY <- 1e-3
## END TEST ##

####################
## Define options ##
####################


###############
## Load data ##
###############

# Load sample metadata
sample_metadata <- fread(args$metadata) %>% .[pass_rnaQC==T & sample%in%args$samples]
print(table(sample_metadata$class))

# Load SingleCellExperiment
sce <- load_SingleCellExperiment(args$sce, cells = sample_metadata$cell, normalise = T)
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

# Load gene metadata
gene_metadata <- fread(io$gene_metadata) %>% 
  .[!grepl("^[Gm|Rik]",symbol)] %>%
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
# genes.chrY <- genes.chrY[!genes.chrY=="ENSMUSG00000096768"]

dt <- args$samples %>% map(function(i) {
  sce.filt <- sce[,sce$sample==i] %>% .[c(genes.chrY,genes.chrX,genes.chr10),]
  data.table(
    sample = i,
    symbol = rownames(sce.filt),
    expr = logcounts(sce.filt) %>% Matrix::rowMeans()
  )
}) %>% rbindlist %>% merge(gene_metadata[,c("chr","ens_id","symbol")], by="symbol")


##########
## Plot ##
##########

to.plot <- dt[chr=="chrY"] %>% 
 .[,foo:=mean(expr),by="symbol"] %>% .[foo>0] %>% .[,foo:=NULL] 

p <- ggbarplot(to.plot, x="symbol", y="expr", facet="sample", fill="gray70") +
# p <- ggbarplot(to.plot, x="ens_id", y="counts", facet="sample", fill="gray70") +
  labs(x="", y="Expression levels") +
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
  .[,.(expr=mean(expr)),by=c("sample","chr")] %>%
  dcast(sample~chr, value.var="expr") %>%
  .[,ratioY:=chrY/chr10] %>% .[,ratioX:=chrX/chr10] %>%
  .[,sex:=c("female","male")[as.numeric(ratioY>=args$threshold.ratioY)+1]]

p <- ggbarplot(sex_assignment.dt, x="sample", y="ratioY", fill="sex", sort.val = "asc", palette="Dark2") +
  labs(x="", y="chrY/chr1 expr ratio") +
  guides(x = guide_axis(angle = 90)) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(colour="black",size=rel(0.6)),
    axis.text.y = element_text(colour="black",size=rel(0.8))
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
  merge(sex_assignment.dt[,c("sex","sample")], by="sample")

p <- ggbarplot(to.plot, x="sample", y="expr", fill="sex") +
  labs(x="", y="Xist expression") +
  guides(x = guide_axis(angle = 90)) +
  theme(
    axis.text.x = element_text(colour="black",size=rel(0.6)),
    axis.text.y = element_text(colour="black",size=rel(0.8))
  )

pdf(sprintf("%s/pdf/Xist_expr.pdf",args$outdir))
print(p)
dev.off()

###########################################
## Scatterplot of chrY/chr1 vs Xist expr ##
###########################################

to.plot <- dt[symbol=="Xist"] %>% 
  setnames("expr","Xist_expr") %>%
  merge(sex_assignment.dt, by="sample")

ggscatter(to.plot, x="ratioY", y="Xist_expr", shape=21, fill="sex", size=2.5) +
  labs(x="chrX/chr1 ratio", y="Xist expression") +
  theme(
    axis.text = element_text(colour="black",size=rel(0.8)),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )

############################
## Update sample metadata ##
############################

sample_metadata_after_sex_assignment <- fread(io$metadata) %>% 
  merge(sex_assignment.dt[,c("sample","sex")], by="sample", all.x = TRUE)

##########
## Save ##
##########

fwrite(sample_metadata_after_sex_assignment, file.path(args$outdir,"sample_metadata_after_sex_assignment.txt.gz"), sep="\t", quote=F, na="NA")
fwrite(sex_assignment.dt, file.path(args$outdir,"sex_assignment.txt.gz"))





