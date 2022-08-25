
source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$outdir <- file.path(io$basedir,"results/individual_genes/pseudobulk/per_sample/imprinting"); dir.create(io$outdir, showWarnings = F)
io$sce.pseudobulk <- file.path(io$basedir,"results/pseudobulk/SingleCellExperiment_pseudobulk_class_sample_celltype.rds")

# Define options
opts$celltypes = c(
  "Spinal_cord",
  "Blood_progenitors_2",
  "Cardiomyocytes",
  "NMP",
  "ExE_endoderm",
  # "ExE_ectoderm"
  "Visceral_endoderm"
)

opts$classes <- c(
  "WT",
  "Dnmt3a_KO",
  "Dnmt3b_KO",
  "Dnmt1_KO"
)

opts$min.cells <- 10

#####################
## Load imprinting ##
#####################

imprinting.dt <- fread(io$imprinted.genes) %>%
  setnames(c("gene","allele"))

#########################
## Load RNA expression ##
#########################

# Load SingleCellExperiment object
sce <- readRDS(io$sce.pseudobulk)

# Add metadata
sce$class <- stringr::str_split(colnames(sce), pattern = "-") %>% map_chr(1)
sce$sample <- stringr::str_split(colnames(sce), pattern = "-") %>% map_chr(2)
sce$celltype <- stringr::str_split(colnames(sce), pattern = "-") %>% map_chr(3)

# Filter by minimum number of cells
metadata(sce)$n_cells <- metadata(sce)$n_cells[metadata(sce)$n_cells>=opts$min.cells]
sce <- sce[,names(metadata(sce)$n_cells)]

# Filter classes and celltypes manually
sce <- sce[,sce$celltype%in%opts$celltypes]
sce <- sce[,sce$class%in%opts$classes]
metadata(sce)$n_cells <- metadata(sce)$n_cells[stringr::str_split(names(metadata(sce)$n_cells), pattern = "-") %>% map_chr(3) %in% opts$celltypes]
metadata(sce)$n_cells <- metadata(sce)$n_cells[stringr::str_split(names(metadata(sce)$n_cells), pattern = "-") %>% map_chr(1) %in% opts$classes]

# Filter celltypes with not enough observations
stats.dt <- data.table(class=sce$class, celltype=sce$celltype, sample=sce$sample, id=colnames(sce))
sce <- sce[,stats.dt$id]

table(stats.dt$class,stats.dt$celltype)

#######################
## Plotting function ##
#######################

plotting_fn <- function(sce, gene, classes, celltypes, max.expr=NULL, add_number_observations=TRUE, scale_expr=TRUE) {
  
  sce.filt <- sce[gene,sce$celltype%in%celltypes & sce$class%in%classes]
  
  to.plot <- data.table(
    expr = logcounts(sce.filt)[gene,],
    class = sce.filt$class,
    celltype = sce.filt$celltype,
    sample = sce.filt$sample
  ) %>% .[,celltype:=factor(celltype,levels=celltypes)] %>% 
    .[,class:=factor(class,levels=opts$classes)]
  
  if (!is.null(max.expr)) to.plot[expr>=max.expr,expr:=max.expr]
  to.plot[expr==0,expr:=0.10]
  
  if (scale_expr) { to.plot[,expr:=expr/max(expr)] }
  to.plot.means <- to.plot[,.(expr=mean(expr),sd=sd(expr)), by=c("class","celltype")]
  
  p <- ggplot(to.plot.means, aes(x=class, y=expr, fill=class)) +
    geom_bar(stat="identity", color="black", width=0.80) +
    geom_jitter(size=1, alpha=0.50, width=0.1, shape=21, data=to.plot) +
    geom_errorbar(aes(ymin=expr-sd, ymax=expr+sd), width=0.25, alpha=0.75, size=0.5) +
    facet_wrap(~celltype, scales="fixed", nrow=2) +
    scale_fill_manual(values=opts$classes.colors) +
    theme_classic() +
    labs(x="",y=sprintf("%s expression",gene)) +
    # guides(x = guide_axis(angle = 90)) +
    theme(
      strip.background = element_blank(),
      plot.title = element_text(hjust = 0.5),
      strip.text = element_text(size=rel(0.85)),
      axis.text.x = element_text(colour="black",size=rel(0.75)),
      axis.text.y = element_text(colour="black",size=rel(0.9)),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(colour="black",size=rel(1.0)),
      legend.position = "none"
      # legend.title = element_blank(),
      # legend.text = element_text(size=rel(0.85))
    )
  
  if (add_number_observations) {
    tmp <- max(to.plot$expr)+0.05
    p <- p + stat_summary(fun.data = function(x){ return(c(y = tmp, label = length(x))) }, geom = "text", size=2.75, data=to.plot)
  }
  
  if (scale_expr) { p <- p + scale_y_continuous(breaks=c(0,1)) }
  
  return(p)
}

################
## Plot genes ##
################

# genes.to.plot <- c("Tex19.1","Rhox6","Rhox9","Trap1a","Xlr3a","Phlda2")
genes.to.plot <- c("Cdkn1c","Grb10","Impact","Peg3","Rhox5","Phlda2")
# genes.to.plot <- imprinting.dt$gene

for (i in genes.to.plot) {
  pdf(sprintf("%s/%s_imprinting_barplot_pseudobulk_per_sample.pdf",io$outdir,i), width=8, height=4)
  p <- plotting_fn(sce, gene = i, celltypes = opts$celltypes, classes = opts$classes) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  print(p)
  dev.off()
}
