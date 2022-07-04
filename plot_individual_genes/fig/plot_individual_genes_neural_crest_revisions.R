here::i_am("plot_individual_genes/plot_individual_genes.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O ##
io$outdir <- file.path(io$basedir,"results/individual_genes/figs"); dir.create(io$outdir, showWarnings = F)

## Define options ##

# Define cell types to plot

# Define classes to plot
opts$classes <- c(
  "WT", 
  # "Dnmt3a_KO", 
  # "Dnmt3b_KO",
  "Dnmt1_KO"
  # "Dnmt3ab_KO"
)

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

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  .[pass_rnaQC==TRUE & celltype.mapped%in%opts$celltypes & class%in%opts$classes] %>%
  .[,celltype:=stringr::str_replace_all(celltype.mapped,opts$rename_celltypes)]

table(sample_metadata$class)
table(sample_metadata$celltype)

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- load_SingleCellExperiment(io$sce, cells=sample_metadata$cell, normalise = TRUE)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

#######################
## Plotting function ##
#######################

plotting_fn <- function(sce, gene, classes, celltypes, max.expr=NULL) {
  
  sce.filt <- sce[gene,sce$celltype%in%celltypes & sce$class%in%classes]
  
  to.plot <- data.table(
    cell = colnames(sce.filt),
    expr = logcounts(sce.filt)[1,],
    class = factor(sce.filt$class, levels=classes),
    celltype = factor(sce.filt$celltype, levels=celltypes)
  )
  
  # For viz purposes
  if (!is.null(max.expr)) {
    to.plot[expr>=max.expr,expr:=max.expr]
  }
  to.plot.jitter <- to.plot[,.SD[sample.int(100)],by=c("class","celltype")]
  
  ggplot(to.plot, aes(x=class, y=expr, fill=class)) +
    geom_violin(scale = "width", alpha=0.8) +
    geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.8) +
    # stat_summary(fun.data = give.n, geom = "text", size=3) +
    geom_jitter(size=1, shape=21, stroke=0.15, alpha=0.4, data=to.plot.jitter, width=0.1) +
    scale_fill_manual(values=opts$classes.colors) +
    facet_wrap(~celltype, scales="fixed", nrow=1) +
    theme_classic() +
    labs(x="",y=sprintf("%s expression",gene)) +
    theme(
      strip.text = element_text(size=rel(0.85)),
      strip.background = element_blank(),
      axis.text.x = element_text(colour="black",size=rel(0.75)),
      axis.text.y = element_text(colour="black",size=rel(0.8)),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(colour="black",size=rel(1.0)),
      legend.position = "none"
    )
}

#############################
## Plot one gene at a time ##
#############################

genes.to.plot <- c("Tex19.1","Rhox6","Rhox9","Trap1a","Xlr3a","Phlda2")
genes.to.plot <- c("Foxd3","Dlx2","Sox10","Sox9","Tfap2c","Tfap2a","Tfap2a")
celltypes.to.plot <- c("Neural_crest")
# celltypes.to.plot <- c( "Gut", "Cardiomyocytes", "Paraxial_mesoderm", "Spinal_cord")

for (i in genes.to.plot) {
  pdf(sprintf("%s/Neural_crest_%s.pdf",io$outdir,i), width=6, height=4)
  print(plotting_fn(sce, gene = i, celltypes = celltypes.to.plot, classes = opts$classes))
  dev.off()
}


##################################
## Plot multiple genes together ##
##################################