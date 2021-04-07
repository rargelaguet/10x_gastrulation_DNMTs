suppressMessages(library(scater))
suppressMessages(library(edgeR))
suppressMessages(library(argparse))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--groupA',    type="character",    help='group A ("class" column in metadata)')
p$add_argument('--groupB',    type="character",    help='group B ("class" column in metadata)')
p$add_argument('--celltype',  type="character",    nargs="+", help='Cell type')
p$add_argument('--test_mode', action="store_true", help='Test mode? subset number of cells')
p$add_argument('--outfile',   type="character",    help='Output file')
args <- p$parse_args(commandArgs(TRUE))

#########
## I/O ##
#########

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
  source("/Users/ricard/10x_gastrulation_DNMTs/utils.R")
  source("/Users/ricard/10x_gastrulation_DNMTs/differential/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/10x_gastrulation_DNMTs/settings.R")
  source("/homes/ricard/10x_gastrulation_DNMTs/utils.R")
  source("/homes/ricard/10x_gastrulation_DNMTs/differential/utils.R")
} else {
  stop("Computer not recognised")
}

# START TEST
# args$groupA <- c("E8.5_Dnmt3aKO_Dnmt3bWT")
# args$groupB <- c("E8.5_WT")
# args$celltype <- "ExE_endoderm"
# args$test_mode <- FALSE
# END TEST

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

############################
## Update sample metadata ##
############################

sample_metadata <- fread(io$metadata) %>%
  .[pass_QC==TRUE & class%in%opts$groups & celltype.mapped%in%args$celltype] %>%
  setnames("class","group") %>% #.[,c("cell","group")] %>%
  .[,group:=factor(group,levels=opts$groups)] %>% 
  setorder(group)


if (isTRUE(args$test_mode)) {
  print("Testing mode activated")
  sample_metadata <- sample_metadata %>% split(.,.$group) %>% map(~ head(.,n=250)) %>% rbindlist
}
table(sample_metadata$group)

if (any(table(sample_metadata$group)<opts$min.cells)) {
  stop("Not enough cells for differential testing")
}

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- load_SingleCellExperiment(io$sce, cells = sample_metadata$cell, normalise = TRUE, remove_non_expressed_genes = FALSE)
sce$group <- sample_metadata$group

# Load gene metadata
gene_metadata <- fread(io$gene_metadata) %>%
  .[symbol%in%rownames(sce)] %>%
  .[,c("symbol","ens_id")] %>% 
  setnames("symbol","gene")

################
## Parse data ##
################

# calculate detection rate per gene
cdr.dt <- data.table(
  gene = rownames(sce),
  detection_rate_groupA = rowMeans(logcounts(sce[,sce$group==opts$groups[1]])>0),
  detection_rate_groupB = rowMeans(logcounts(sce[,sce$group==opts$groups[2]])>0)
)# %>% setnames(c("ens_id",sprintf("detection_rate_%s",opts$groups[1]),sprintf("detection_rate_%s",opts$groups[2])))

# Filter genes
sce <- sce[rownames(sce)%in%gene_metadata$gene,]

################################################
## Differential expression testing with edgeR ##
################################################

out <- doDiffExpr(sce, opts$groups, opts$min_detection_rate_per_group) %>%
  # Add sample statistics
  .[,c("groupA_N","groupB_N"):=list(table(sample_metadata$group)[1],table(sample_metadata$group)[2])]%>% 
  # setnames(c("groupA_N","groupB_N"),c(sprintf("N_%s",opts$groups[1]),sprintf("N_%s",opts$groups[2]))) %>%
  # Add gene statistics
  merge(cdr.dt, all.y=T, by="gene") %>%
  merge(gene_metadata, all.y=T, by="gene") %>%
  # Calculate statistical significance
  .[, sig := (padj_fdr<=opts$threshold_fdr & abs(logFC)>=opts$min.logFC)] %>%
  .[is.na(sig),sig:=FALSE] %>%
  setorder(-sig, padj_fdr, na.last=T)

# Parse columns
out[,c("p.value","padj_fdr","logFC","log_padj_fdr"):=list(signif(p.value,digits=3), signif(padj_fdr,digits=3), round(logFC,3),round(log_padj_fdr,3))]  

##################
## Save results ##
##################

fwrite(out, args$outfile, sep="\t", na="NA", quote=F)
