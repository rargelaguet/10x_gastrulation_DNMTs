here::i_am("plot_individual_genes/pseudobulk/plot_individual_genes_pseudobulk.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$outdir <- file.path(io$basedir,"results_all/individual_genes/pseudobulk"); dir.create(io$outdir, showWarnings = F)
io$sce.pseudobulk <- file.path(io$basedir,"results_all/pseudobulk/SingleCellExperiment_pseudobulk_class_celltype_dataset.rds")

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

# opts$rename_celltypes <- c(
#   "Erythroid3" = "Erythroid",
#   "Erythroid2" = "Erythroid",
#   "Erythroid1" = "Erythroid",
#   "Blood_progenitors_1" = "Blood_progenitors",
#   "Blood_progenitors_2" = "Blood_progenitors",
#   "Anterior_Primitive_Streak" = "Primitive_Streak"
# )

opts$min.cells <- 30

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- readRDS(io$sce.pseudobulk)

# Add metadata
sce$class <- stringr::str_split(colnames(sce), pattern = "-") %>% map_chr(1)
sce$celltype <- stringr::str_split(colnames(sce), pattern = "-") %>% map_chr(2)
sce$dataset <- stringr::str_split(colnames(sce), pattern = "-") %>% map_chr(3)

# Filter by minimum number of cells
sce <- sce[,names(which(metadata(sce)$n_cells>=opts$min.cells))]

# Filter
sce <- sce[,sce$celltype%in%opts$celltypes]
sce <- sce[,sce$class%in%opts$classes]

#############
## Barplots ##
#############

# genes.to.plot <- c("Tex19.1","Morc1","Dppa3","Rex1","Dppa5a","Dppa4","Dppa2","Zfp981")
genes.to.plot <- c("Trap1a","Tex19.1","Rhox5")
# genes.to.plot <- grep("^Hox",rownames(sce),value=T)
# genes.to.plot <- fread(io$atlas.marker_genes)[,gene] %>% unique %>% head(n=3)

for (i in 1:length(genes.to.plot)) {
  
  gene <- genes.to.plot[i]
  
  if (gene %in% rownames(sce)) {
    print(sprintf("%s/%s: %s",i,length(genes.to.plot),gene))
    outfile <- sprintf("%s/%s_barplot_pseudobulk.pdf",io$outdir,gene)
  
    to.plot <- data.table(
      sample = colnames(sce),
      expr = logcounts(sce)[gene,],
      class = sce$class,
      celltype = sce$celltype,
      dataset = sce$dataset
    ) %>% .[,dataset:=factor(dataset,levels=c("KO","CRISPR"))] %>% .[,celltype:=factor(celltype,levels=opts$celltypes)] %>% .[,class:=factor(class,levels=opts$classes)]
      
    p_list <- list()
    for (j in levels(to.plot$dataset)) {
      
      p_list[[j]] <- ggplot(to.plot[dataset==j], aes(x=class, y=expr, fill=class)) +
        geom_bar(stat="identity", color="black", width=0.75) +
        facet_wrap(~celltype, scales="fixed") +
        # scale_fill_brewer(palette="Dark2") +
        scale_fill_manual(values=opts$classes.colors) +
        theme_classic() +
        labs(x="",y=sprintf("%s expression",gene), title=j) +
        guides(x = guide_axis(angle = 90)) +
        theme(
          plot.title = element_text(hjust = 0.5),
          strip.text = element_text(size=rel(0.85)),
          axis.text.x = element_text(colour="black",size=rel(0.9)),
          axis.text.y = element_text(colour="black",size=rel(0.9)),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(colour="black",size=rel(1.0)),
          legend.position = "none"
          # legend.title = element_blank(),
          # legend.text = element_text(size=rel(0.85))
        )
    }
    pdf(outfile, width=17, height=9)
    print(cowplot::plot_grid(plotlist=p_list, ncol = 2))
    dev.off()
        
  } else {
    print(sprintf("%s not found",gene))
  }
}

