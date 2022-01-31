here::i_am("plot_individual_genes/pseudobulk/plot_individual_genes_pseudobulk_paga.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$outdir <- file.path(io$basedir,"results/individual_genes/pseudobulk/paga"); dir.create(io$outdir, showWarnings = F)
io$sce.pseudobulk <- file.path(io$basedir,"results/pseudobulk/SingleCellExperiment_pseudobulk_class_celltype.rds")
# io$sce.pseudobulk <- file.path(io$basedir,"results/pseudobulk/SingleCellExperiment_pseudobulk_class_celltype_dataset.rds")

# Define options
opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  # "PGC",
  "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  "Mixed_mesoderm",
  "Intermediate_mesoderm",
  "Caudal_Mesoderm",
  "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "Pharyngeal_mesoderm",
  "Cardiomyocytes",
  "Allantois",
  "ExE_mesoderm",
  "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  # "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm"
  # "ExE_ectoderm",
  # "Parietal_endoderm"
)

opts$classes <- c(
  # "Dnmt3a_HET_Dnmt3b_KO",
  # "Dnmt3a_HET_Dnmt3b_WT",
  # "Dnmt3a_KO_Dnmt3b_HET",
  "WT",
  "Dnmt3a_KO",
  "Dnmt3b_KO",
  "Dnmt1_KO",
  "Dnmt3ab_KO"
)

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- readRDS(io$sce.pseudobulk)

# Add metadata
sce$class <- stringr::str_split(colnames(sce), pattern = "-") %>% map_chr(1)
sce$celltype <- stringr::str_split(colnames(sce), pattern = "-") %>% map_chr(2)
# sce$dataset <- stringr::str_split(colnames(sce), pattern = "-") %>% map_chr(3)

# Filter by minimum number of cells
# sce <- sce[,names(which(metadata(sce)$n_cells>=opts$min.cells))]

# Filter
sce <- sce[,sce$celltype%in%opts$celltypes]
sce <- sce[,sce$class%in%opts$classes]

opts$celltypes <- unique(sce$celltype)
###############
## Load PAGA ##
###############

source(here::here("load_paga_graph.R"))

# Plot graph structure
p <- ggnet2(
  net = net.paga,
  mode = c("x", "y"),
  node.size = 0,
  edge.size = 0.15,
  edge.color = "grey",
  label = FALSE,
  label.size = 2.3
)

paga.celltype.order <- p$data$label

##########
## Plot ##
##########

# Define color scale
rna.col.seq <- chromvar.col.seq <- round(seq(0,1,0.1), 2)
rna.colors <- colorRampPalette(c("gray92", "darkgreen"))(length(rna.col.seq))

# Define genes to plot
# genes.to.plot <- rownames(rna.sce)[grep("Gata",rownames(rna.sce))]
genes.to.plot <- c("Foxa2","Tfap2a","Mesp1")

# i <- "Foxa2"; j <- "Dnmt1_KO"
for (i in genes.to.plot) {
  print(i)
  
  max.expr <- max(logcounts(sce[i,]))
  p_list <- list()
  for (j in opts$classes) {
    # Subset KO samples from the SingleCellExperiment
    sce.tmp <- sce[,sce$class==j]
    metadata(sce.tmp)$n_cells <- metadata(sce.tmp)$n_cells[stringr::str_split(names(metadata(sce.tmp)$n_cells), pattern = "-") %>% map_chr(1) == j]
    names(metadata(sce.tmp)$n_cells) <- stringr::str_split(names(metadata(sce.tmp)$n_cells), pattern = "-") %>% map_chr(2)
    colnames(sce.tmp) <- sce.tmp$celltype
    
    expr.values <- rep(as.numeric(NA),length(paga.celltype.order)); names(expr.values) <- paga.celltype.order
    expr.values[colnames(sce.tmp)] <- logcounts(sce.tmp[i,])[1,] / max.expr 
    # expr.values <- logcounts(sce.tmp[i,])[1,] %>% minmax.normalisation()
    
    expr.colors <- round(expr.values,1) %>% map(~ rna.colors[which(rna.col.seq == .)]) %>% unlist
    
    p_list[[j]] <- p + geom_text(label = "\u25D0", aes(x=x, y=y), color=expr.colors, size=15, family = "Arial Unicode MS",
                           data = p$data[p$data$label%in%names(expr.colors),c("x","y")] %>% dplyr::mutate(expr=expr.colors)) +
      scale_colour_manual(values=expr.colors) + 
      labs(title=j) +
      theme(
        plot.title = element_text(hjust = 0.5)
      )
    
  }
  png(sprintf("%s/%s_pseudobulk_paga.png",io$outdir,i), width = 1000, height = 300)
  print(cowplot::plot_grid(plotlist=p_list, nrow=1))
  dev.off()
  # # pdf(sprintf("%s/%s_rna_expression_paga.pdf",io$outdir,i), width=5, height=3.5)
}
