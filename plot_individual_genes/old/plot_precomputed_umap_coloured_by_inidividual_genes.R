suppressPackageStartupMessages(library(SingleCellExperiment))

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation10x/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation10x/settings.R")
}
io$outdir <- paste0(io$basedir,"/results/individual_genes/umap")

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- readRDS(io$rna.sce)[,as.character(sample_metadata$cell)]

# Remove genes that are not expressed
sce <- sce[rowMeans(counts(sce))>0,]

# Load gene metadata
gene_metadata <- fread(io$gene_metadata) %>%
  .[ens_id%in%rownames(sce)]

################
## Parse data ##
################

# Rename genes
foo <- gene_metadata$symbol
names(foo) <- gene_metadata$ens_id
sce <- sce[rownames(sce) %in% names(foo),]
rownames(sce) <- foo[rownames(sce)]


# Extract UMAP coordinates from the metadata
umap.dt <- sample_metadata[,c("umapX","umapY")]

###################################
## Plot dimensionality reduction ##
###################################

genes.to.plot <- c("Dnmt3l","Apoe")

for (i in 1:length(genes.to.plot)) {
  gene <- genes.to.plot[i]
  print(sprintf("%s/%s: %s",i,length(genes.to.plot),gene))
  
  # Create data.table to plot
  to.plot <- data.table(
    cell = colnames(sce),
    expr = logcounts(sce[i,])[1,]
  ) %>% cbind(umap.dt)

  # Plot
  p <- ggplot(to.plot, aes(x=umapX, y=umapY, color=expr)) +
    scale_color_gradient(low = "gray80", high = "red") +
    geom_point(size=0.25) +
    theme_classic() +
    theme(
      legend.title = element_blank(),
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank()
    )
    
  # pdf(sprintf("%s/%s.pdf",io$outdir,i), width=8, height=3.5, useDingbats = F)
  # ggsave("ggtest.png", width = 3.25, height = 3.25, dpi = 1200)
  jpeg(sprintf("%s/%s_umap.jpeg",io$outdir,gene), width = 600, height = 600)
  print(p)
  dev.off()
}
