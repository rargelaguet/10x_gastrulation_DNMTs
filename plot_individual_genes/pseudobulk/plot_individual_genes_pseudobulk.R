here::i_am("plot_individual_genes/pseudobulk/plot_individual_genes_pseudobulk.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$outdir <- file.path(io$basedir,"results_new/individual_genes/pseudobulk")
io$sce.pseudobulk <- file.path(io$basedir,"results_new/pseudobulk/SingleCellExperiment_pseudobulk_class_celltype.rds")

# Define options
opts$min.cells <- 25

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

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- readRDS(io$sce.pseudobulk)

# Add metadata
sce$class <- stringr::str_split(colnames(sce), pattern = "-") %>% map_chr(1)
sce$celltype <- stringr::str_split(colnames(sce), pattern = "-") %>% map_chr(2)

# Filter celltypes
sce <- sce[,sce$celltype%in%opts$celltypes]

# Filter by minimum number of cells
sce <- sce[,names(which(metadata(sce)$n_cells>=opts$min.cells))]

#############
## Barplots ##
#############

genes.to.plot <- c("Tex19.1","Morc1","Dppa3","Rex1","Dppa5a","Dppa4","Dppa2","Zfp981")
genes.to.plot <- c("Pou5f1","Epcam","Fst","Lefty2","Cdkn1c","Acta1")
genes.to.plot <- grep("^Hox",rownames(sce),value=T)

for (i in 1:length(genes.to.plot)) {
  
  gene <- genes.to.plot[i]
  
  if (gene %in% rownames(sce)) {
    print(sprintf("%s/%s: %s",i,length(genes.to.plot),gene))
    outfile <- sprintf("%s/%s.pdf",io$outdir,gene)
  
    to.plot <- data.table(
      sample = colnames(sce),
      expr = logcounts(sce)[gene,],
      class = sce$class,
      celltype = sce$celltype
    )
    
    p <- ggplot(to.plot, aes(x=class, y=expr, fill=class)) +
      geom_bar(stat="identity", color="black") +
      scale_fill_brewer(palette="Dark2") +
      facet_wrap(~celltype, scales="fixed") +
      theme_classic() +
      labs(x="",y=sprintf("%s expression",gene)) +
      guides(x = guide_axis(angle = 90)) +
      theme(
        strip.text = element_text(size=rel(0.85)),
        axis.text.x = element_text(colour="black",size=rel(0.9)),
        axis.text.y = element_text(colour="black",size=rel(0.9)),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(colour="black",size=rel(1.0)),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=rel(0.85))
      )
    
      pdf(outfile, width=10, height=9)
      print(p)
      dev.off()
        
  } else {
    print(sprintf("%s not found",gene))
  }
}

