here::i_am("plot_individual_genes/pseudobulk/plot_individual_genes_pseudobulk.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$outdir <- file.path(io$basedir,"results/individual_genes/pseudobulk/per_sample"); dir.create(io$outdir, showWarnings = F)
# io$sce.pseudobulk <- file.path(io$basedir,"results/pseudobulk/SingleCellExperiment_pseudobulk_class_celltype_dataset.rds")
io$sce.pseudobulk <- file.path(io$basedir,"results/pseudobulk/SingleCellExperiment_pseudobulk_class_sample_celltype.rds")

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
  "Dnmt1_KO"
  # "Dnmt3ab_KO"
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
stats.dt <- data.table(class=sce$class, celltype=sce$celltype, sample=sce$sample, id=colnames(sce)) %>% 
  .[,N:=length(unique(class)),by=c("celltype")] %>% .[N>=3]
sce <- sce[,stats.dt$id]

#################################
## Barplots one gene at a time ##
#################################

# genes.to.plot <- c("Bhmt", "Spic", "Trim43a", "Pramel6", "Nanog", "Khdc3", "Tfap2c", "Calcoco2", "Dppa3", "Dppa4", "Esrrb", "Crxos")
# genes.to.plot <- c("Bex1", "Cldn6", "Dstn", "Esrrb", "Etv5", "Fbxo15", "Fgf4", "Pou5f1", "Gdf3", "Idh1", "Jarid2", "Klf2", "Klf4", "Klf5", "Lats2", "Lgals1", "Kdm1a", "Myh9", "Nanog", "Nr0b1", "Pdk1", "Peg3", "Pfkp", "Phldb2", "Plk2", "Rif1", "Nacc1", "Zfp281", "Zfp42", "Sall1", "Sall4", "Sox2", "Spp1", "Dppa3", "Tagln", "Tbx3", "Tcl1", "Tdgf1", "Sfrp1", "Tdh", "Tead1", "Tet1", "Tet2", "Tfcp2l1", "Thbs1", "Tuba1a", "Tubb6", "Upp1", "Utf1", "Wbp5", "Foxa2", "Cer1", "T", "Nodal", "Lef1", "Fgf5")
# genes.to.plot <- c("Dnmt3b", "Gng3", "Nr0b1", "Pecam1")
# genes.to.plot <- fread("/Users/argelagr/data/10x_gastrulation_DNMTs/results/differential/Dnmt1_KO/Neural_crest_WT_vs_Dnmt1_KO.txt.gz") %>% .[padj_fdr<0.01 & abs(logFC)>1.5,gene]
genes.to.plot <- grep("^Hox",rownames(sce),value=T)
# genes.to.plot <- fread(io$atlas.marker_genes)[,gene] %>% unique %>% head(n=3)
# genes.to.plot <- fread(io$atlas.marker_genes) %>% .[grep("Neural_crest",celltype),gene] %>% unique

# genes.to.plot <- genes.to.plot[genes.to.plot%in%fread(io$atlas.marker_genes)[,gene]]

genes.to.plot <- genes.to.plot[genes.to.plot%in%rownames(sce)]

for (i in 1:length(genes.to.plot)) {
  
  gene <- genes.to.plot[i]

  print(sprintf("%s/%s: %s",i,length(genes.to.plot),gene))
  outfile <- sprintf("%s/%s_barplot_pseudobulk_per_sample.pdf",io$outdir,gene)
  
  # Filter values with small N
  to.plot <- data.table(
    # id = colnames(sce),
    expr = logcounts(sce)[gene,],
    class = sce$class,
    celltype = sce$celltype,
    sample = sce$sample
  ) %>% 
    .[,celltype:=factor(celltype,levels=opts$celltypes)] %>% 
    .[,class:=factor(class,levels=opts$classes)]
    
  to.plot[expr==0,expr:=0.10]
  to.plot.means <- to.plot[,.(expr=mean(expr)), by=c("class","celltype")]
  
  p <- ggplot(to.plot.means, aes(x=class, y=expr, fill=class)) +
    geom_bar(stat="identity", color="black", width=0.65) +
    geom_jitter(size=1.25, alpha=0.75, width=0.15, shape=21, data=to.plot) +
    facet_wrap(~celltype, scales="fixed") +
    # scale_fill_brewer(palette="Dark2") +
    scale_fill_manual(values=opts$classes.colors) +
    theme_classic() +
    labs(x="",y=sprintf("%s expression",gene)) +
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
  pdf(outfile, width=14, height=9)
  print(p)
  dev.off()
}


#######################################
## Barplots multiple genes at a time ##
#######################################

# genes.to.plot <- c("Trap1a","Tex19.1","Rhox5")
# 
# p <- ggplot(to.plot[dataset==j], aes(x=class, y=expr, fill=class)) +
#   geom_bar(stat="identity", color="black", width=0.75) +
#   facet_wrap(~celltype, scales="fixed") +
#   # scale_fill_brewer(palette="Dark2") +
#   scale_fill_manual(values=opts$classes.colors) +
#   theme_classic() +
#   labs(x="",y=sprintf("%s expression",gene), title=j) +
#   guides(x = guide_axis(angle = 90)) +
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     strip.text = element_text(size=rel(0.85)),
#     axis.text.x = element_text(colour="black",size=rel(0.9)),
#     axis.text.y = element_text(colour="black",size=rel(0.9)),
#     axis.ticks.x = element_blank(),
#     axis.title.y = element_text(colour="black",size=rel(1.0)),
#     legend.position = "none"
#     # legend.title = element_blank(),
#     # legend.text = element_text(size=rel(0.85))
#   )
# 
# pdf(file.path(io$outdir,"logFC_heatmap_pseudobulk_ExE_genes.pdf"), width=8, height=5)
# print(p)
# dev.off()
