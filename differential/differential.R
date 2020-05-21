######################################################################################
## Script to do differential expression between WT and KO, separately per cell type ##
######################################################################################

suppressMessages(library(SingleCellExperiment))
suppressMessages(library(scater))
suppressMessages(library(edgeR))
suppressMessages(library(argparse))

## Initialize argument parser ##
p <- ArgumentParser(description='')
p$add_argument('--groupA',    type="character",    help='group A')
p$add_argument('--groupB',    type="character",    help='group B')
p$add_argument('--test',      type="character",    help='Statistical test')
p$add_argument('--celltype',  type="character",    nargs="+", help='Cell type')
p$add_argument('--test_mode', action="store_true", help='Test mode? subset number of cells')
p$add_argument('--outfile',   type="character",    help='Output file')
args <- p$parse_args(commandArgs(TRUE))

# Sanity checks
stopifnot(args$test%in%c("edgeR","t-test","wilcoxon"))

#########
## I/O ##
#########

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
  source("/Users/ricard/10x_gastrulation_DNMTs/differential/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/10x_gastrulation_DNMTs/settings.R")
  source("/homes/ricard/10x_gastrulation_DNMTs/differential/utils.R")
} else {
  stop("Computer not recognised")
}

################
## START TEST ##
################

# args$groupA <- c("Dnmt3aWT_Dnmt3bHET")
# args$groupB <- c("Dnmt3aWT_Dnmt3bKO")
# args$celltype <- opts$celltypes
# args$test <- c("edgeR")
# args$test_mode <- TRUE

##############
## END TEST ##
##############

#############
## Options ##
#############

# Define cell types
opts$celltype <- args$celltype
# stopifnot(opts$celltype %in% opts$celltype.1)

# Define groups
opts$groups <- c(args$groupA,args$groupB)

# Define FDR threshold
opts$threshold_fdr <- 0.01

# Define minimum logFC for significance
opts$min.logFC <- 1

# Minimum number of cells per group for differential testing
opts$min.cells <- 25

# For a given gene, the minimum fraction of cells that must express it in at least one group
opts$min_detection_rate_per_group <- 0.40

###############
## Load data ##
###############

stopifnot(all(opts$groups%in%unique(sample_metadata$class)))

# Update cell metadata
sample_metadata <- sample_metadata %>%
  .[class%in%opts$groups & celltype.mapped%in%args$celltype] %>%
  setnames("class","group") %>%
  .[,c("cell","group")]

# Sort cells so that groupA comes before groupB
sample_metadata[,group:=factor(group,levels=opts$groups)] %>% setorder(group)

if (isTRUE(args$test_mode)) {
  print("Testing mode activated")
  sample_metadata <- sample_metadata %>% split(.,.$group) %>% map(~ head(.,n=250)) %>% rbindlist
}
table(sample_metadata$group)

if (any(table(sample_metadata$group)<opts$min.cells)) {
  stop("Not enough cells for differential testing")
}

# Load SingleCellExperiment object
sce <- readRDS(io$sce)[,sample_metadata$cell]
sce$group <- sample_metadata$group

# Load gene metadata
gene_metadata <- fread(io$gene_metadata) %>%
  .[ens_id%in%rownames(sce)] %>%
  .[,c("symbol","ens_id")] %>% 
  setnames("symbol","gene")

################
## Parse data ##
################

# calculate detection rate per gene
cdr.dt <- data.table(
  ens_id = rownames(sce),
  detection_rate_A = rowMeans(logcounts(sce[,sce$group==opts$groups[1]])>0),
  detection_rate_B = rowMeans(logcounts(sce[,sce$group==opts$groups[2]])>0)
) %>% setnames(c("ens_id",sprintf("detection_rate_%s",opts$groups[1]),sprintf("detection_rate_%s",opts$groups[2])))

# Filter genes
sce <- sce[rownames(sce)%in%gene_metadata$ens_id,]

################################################
## Differential expression testing with edgeR ##
################################################

out <- doDiffExpr(sce, opts$groups, args$test, opts$min_detection_rate_per_group) %>%
  # Add sample statistics
  .[,c("groupA_N","groupB_N"):=list(table(sample_metadata$group)[1],table(sample_metadata$group)[2])]%>% 
  setnames(c("groupA_N","groupB_N"),c(sprintf("N_%s",opts$groups[1]),sprintf("N_%s",opts$groups[2]))) %>%
  # Add gene statistics
  merge(cdr.dt, all.y=T, by="ens_id") %>%
  merge(gene_metadata, all.y=T, by="ens_id") %>%
  # Calculate statistical significance
 .[, sig := (padj_fdr<=opts$threshold_fdr & abs(logFC)>=opts$min.logFC)] %>%
  setorder(-sig, padj_fdr, na.last=T)
  

##################
## Save results ##
##################

fwrite(out, args$outfile, sep="\t", na="NA", quote=F)
