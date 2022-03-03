here::i_am("differential/differential.R")

suppressMessages(library(scater))
suppressMessages(library(edgeR))

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--groupA',    type="character",    help='group A')
p$add_argument('--groupB',    type="character",    help='group B')
p$add_argument('--celltypes',    type="character",    default="all", nargs="+", help='Celltypes to use')
p$add_argument('--group_label',    type="character",    help='Group label')
p$add_argument('--outfile',   type="character",    help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
# args$groupA <- "WT"
# args$groupB <- "Dnmt1_KO"
# args$group_label <- "class"
# args$celltypes <- c("Blood_progenitors")
## END TEST

#####################
## Define settings ##
#####################

# Load utils
source(here::here("differential/utils.R"))

# Sanity checks
# if (args$celltypes[1]=="all") {
#   args$celltypes <- opts$celltypes
# } else {
#   stopifnot(args$celltypes%in%opts$celltypes)
# }

# Define groups
opts$groups <- c(args$groupA,args$groupB)

# Define FDR threshold
opts$threshold_fdr <- 0.01

# Define minimum logFC for significance
opts$min.logFC <- 1

# For a given gene, the minimum fraction of cells that must express it in at least one group
opts$min_detection_rate_per_group <- 0.25

# Rename celltypes
opts$rename_celltypes <- c(
  "Erythroid3" = "Erythroid",
  "Erythroid2" = "Erythroid",
  "Erythroid1" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors"
  # "Intermediate_mesoderm" = "Mixed_mesoderm",
  # "Paraxial_mesoderm" = "Mixed_mesoderm",
  # "Nascent_mesoderm" = "Mixed_mesoderm",
  # "Pharyngeal_mesoderm" = "Mixed_mesoderm"
  # "Visceral_endoderm" = "ExE_endoderm"
)

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,opts$rename_celltypes)] %>%
  .[,dataset:=ifelse(grepl("Grosswendt",sample),"CRISPR","KO")] %>%
  .[pass_rnaQC==TRUE & celltype.mapped%in%args$celltypes]

stopifnot(args$group_label%in%colnames(sample_metadata))

sample_metadata <- sample_metadata %>%
  setnames(args$group_label,"group") %>%
  # .[,group:=eval(as.name(args$group_label))] %>%
  .[group%in%c(args$groupA,args$groupB)] %>%
  .[,group:=factor(group,levels=opts$groups)] %>% setorder(group) # Sort cells so that groupA comes before groupB

table(sample_metadata$group)

#########################
## Load RNA expression ##
#########################

# Load SingleCellExperiment object
sce <- load_SingleCellExperiment(
  file = io$sce, 
  normalise = TRUE, 
  cells = sample_metadata$cell
)
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

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
  rownames(sce),
  rowMeans(logcounts(sce[,sce$group==opts$groups[1]])>0) %>% round(2),
  rowMeans(logcounts(sce[,sce$group==opts$groups[2]])>0) %>% round(2)
) %>% setnames(c("gene",sprintf("detection_rate_%s",opts$groups[1]),sprintf("detection_rate_%s",opts$groups[2])))
# .[,cdr_diff:=abs(out[,(sprintf("detection_rate_%s",opts$groups[1])),with=F][[1]] - out[,(sprintf("detection_rate_%s",opts$groups[2])),with=F][[1]])] %>%

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
  # .[, sig := (padj_fdr<=opts$threshold_fdr & abs(logFC)>=opts$min.logFC)] %>%
  # .[is.na(sig),sig:=FALSE] %>%
  setorder(padj_fdr, na.last=T)

# Parse columns
out[,c("p.value","padj_fdr","logFC","log_padj_fdr"):=list(signif(p.value,digits=3), signif(padj_fdr,digits=3), round(logFC,3),round(log_padj_fdr,3))]

##################
## Save results ##
##################

fwrite(out, args$outfile, sep="\t", na="NA", quote=F)
