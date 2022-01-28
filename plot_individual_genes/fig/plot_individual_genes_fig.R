here::i_am("plot_individual_genes/plot_individual_genes.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O ##
io$outdir <- file.path(io$basedir,"results_all/individual_genes/figs"); dir.create(io$outdir, showWarnings = F)

## Define options ##

# Define cell types to plot

# Define classes to plot
opts$classes <- c(
  "WT", 
  "Dnmt3a_KO", 
  "Dnmt3b_KO",
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
  .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,opts$rename_celltypes)]

table(sample_metadata$class)
table(sample_metadata$celltype.mapped)

# foo <- fread(io$metadata) %>% .[pass_rnaQC==TRUE & celltype.mapped=="Visceral_endoderm"] %>% .[,.N,by="class"]
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
  
  sce.filt <- sce[gene,sce$celltype.mapped%in%celltypes & sce$class%in%classes]
  
  to.plot <- data.table(
    cell = colnames(sce.filt),
    expr = logcounts(sce.filt)[1,],
    class = factor(sce.filt$class, levels=classes),
    celltype = factor(sce.filt$celltype.mapped, levels=celltypes)
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

####################
## Plot ExE genes ##
####################

genes.to.plot <- c("Rhox5","Rhox6","Rhox9","Trap1a","Xlr3a")
celltypes.to.plot <- c("Blood_progenitors","Gut","ExE_ectoderm", "ExE_endoderm")
# celltypes.to.plot <- c( "Gut", "Cardiomyocytes", "Paraxial_mesoderm", "Spinal_cord")

for (gene in genes.to.plot) {
  pdf(sprintf("%s/ExE_%s.pdf",io$outdir,gene), width=8, height=4)
  plotting_fn(sce, gene = gene, celltypes = celltypes.to.plot, classes = opts$classes)
  dev.off()
}

###############################
## Primed Pluripotency genes ##
###############################

# For the primed: Pim2 and Fgf5 are good choices, Gng3 and Slc7a3 may not be known as pluripotency markers. Other good genes to include Lef1, Pou3f1, Sox3, Otx2 (although it might be expressed in other lineages), Foxd3, Sall2, Tead2, Hes6, Nodal. 
# You can also check Etv4, Etv5, Myc, Sox11, Fzd7, Sfrp2, Cdh2 (they are higher in primed than naive, but may not be specific for pluripotency - as many of these genes actually)

celltypes.to.plot <- c("Rostral_neurectoderm","Somitic_mesoderm","Forebrain_Midbrain_Hindbrain", "Spinal_cord","Mesenchyme","Cardiomyocytes")
genes.to.plot <-  c("Pou5f1", "Utf1", "Pim2", "Slc7a3", "Fgf5","Gng3") # "Dnmt3b"

for (gene in genes.to.plot) {
  pdf(sprintf("%s/plurypotency_%s.pdf",io$outdir,gene), width=8, height=4)
  plotting_fn(sce, gene = gene, celltypes = celltypes.to.plot, classes = opts$classes)
  dev.off()
}

################################
## General Pluripotency genes ##
################################

celltypes.to.plot <- c("Rostral_neurectoderm","Somitic_mesoderm","Forebrain_Midbrain_Hindbrain", "Spinal_cord","Mesenchyme","Cardiomyocytes")
genes.to.plot <-  c("Pou5f1", "Sall4", "Utf1", "Tdgf1", "Lin28a", "Zic2", "Zic3", "Nanog", "Sox2")

for (gene in genes.to.plot) {
  pdf(sprintf("%s/plurypotency_%s.pdf",io$outdir,gene), width=8, height=4)
  plotting_fn(sce, gene = gene, celltypes = celltypes.to.plot, classes = opts$classes)
  dev.off()
}

##############################
## Naive Pluripotency genes ##
##############################

celltypes.to.plot <- c("Rostral_neurectoderm","Somitic_mesoderm","Forebrain_Midbrain_Hindbrain", "Spinal_cord","Mesenchyme","Cardiomyocytes")
genes.to.plot <-  c("Esrrb", "Klf2", "Klf4", "Klf5", "Tfcp2l1", "Nr5a2", "Fgf4", "Nr0b1")
# genes.to.plot <-  c("Zfp42", "Dppa5a", "Dppa3", "Dppa4", "Tfap2c", "Pecam1")

for (gene in genes.to.plot) {
  pdf(sprintf("%s/plurypotency_%s.pdf",io$outdir,gene), width=8, height=4)
  plotting_fn(sce, gene = gene, celltypes = celltypes.to.plot, classes = opts$classes)
  dev.off()
}

####################
## Plot Hox genes ##
####################

genes.to.plot <- c("Hoxc9","Hoxc8","Hoxb9","Hoxa9")
celltypes.to.plot <- c("ExE_mesoderm", "Caudal_Mesoderm", "NMP", "Somitic_mesoderm")

for (gene in genes.to.plot) {
  pdf(sprintf("%s/hox_%s.pdf",io$outdir,gene), width=8, height=4)
  plotting_fn(sce, gene = gene, celltypes = celltypes.to.plot, classes = opts$classes)
  dev.off()
}

