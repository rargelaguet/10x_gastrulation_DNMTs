here::i_am("plot_individual_genes/celltype_markers/plot_celltype_markers.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O 
io$outdir <- file.path(io$basedir,"results_new/individual_genes/celltype_markers"); dir.create(io$outdir, showWarnings = F)

opts$classes <- c(
  "E8.5_WT",
  "E8.5_Dnmt3aHET_Dnmt3bWT",
  "E8.5_Dnmt3aKO_Dnmt3bWT",
  "E8.5_Dnmt3aKO_Dnmt3bHET",
  "E8.5_Dnmt3aWT_Dnmt3bKO",
  "E8.5_Dnmt3aHET_Dnmt3bKO",
  "E8.5_Dnmt3aKO_Dnmt3bKO",
  "E8.5_Dnmt1KO"
)

# opts$rename_celltypes <- c(
#   "Erythroid3" = "Erythroid",
#   "Erythroid2" = "Erythroid",
#   "Erythroid1" = "Erythroid",
#   "Blood_progenitors_1" = "Blood_progenitors",
#   "Blood_progenitors_2" = "Blood_progenitors",
#   "Anterior_Primitive_Streak" = "Primitive_Streak"
# )

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  # .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,opts$rename_celltypes)] %>%
  .[pass_rnaQC==TRUE & celltype.mapped%in%opts$celltypes & class%in%opts$classes]

# Only consider cell types with sufficient observations in WT cells
celltypes.to.use <- sample_metadata[class=="E8.5_WT",.(N=.N),by="celltype.mapped"] %>% .[N>=50,celltype.mapped]
sample_metadata <- sample_metadata[celltype.mapped%in%celltypes.to.use]

table(sample_metadata$class)
table(sample_metadata$celltype.mapped)

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- load_SingleCellExperiment(io$sce, cells=sample_metadata$cell, normalise = TRUE)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

#######################
## Load marker genes ##
#######################

marker_genes.dt <- fread(io$atlas.marker_genes) %>%
  .[gene%in%rownames(sce) & celltype%in%celltypes.to.use]

##########
## Plot ##
##########

# i <- "Neural_crest"
for (i in opts$celltypes) {
  genes.to.plot <- marker_genes.dt[celltype==i] %>% setorder(-score) %>% head(n=10) %>% .$gene
  for (j in genes.to.plot) {
    
    outfile <- sprintf("%s/%s_%s.pdf",io$outdir,i,j)
    
    to.plot <- data.table(
      cell = colnames(sce),
      expr = logcounts(sce)[j,]
    ) %>% merge(sample_metadata[,c("cell","sample","class","celltype.mapped")], by="cell") %>%
      .[,N:=.N,by=c("sample","celltype.mapped")] %>% .[N>=15]
    
    p <- ggplot(to.plot, aes(x=class, y=expr, fill=class)) +
      geom_violin(scale = "width", alpha=0.8) +
      geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.8) +
      # stat_summary(fun.data = give.n, geom = "text", size=3) +
      # scale_fill_manual(values=opts$colors) +
      scale_fill_brewer(palette="Dark2") +
      facet_wrap(~celltype.mapped, scales="fixed") +
      theme_classic() +
      labs(x="",y=sprintf("%s expression",j)) +
      guides(x = guide_axis(angle = 90)) +
      theme(
        strip.text = element_text(size=rel(0.85)),
        axis.text = element_text(colour="black",size=rel(0.8)),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(colour="black",size=rel(1.0)),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=rel(0.85))
      )
    
      pdf(outfile, width=10, height=9)
      print(p)
      dev.off()
  }
}
