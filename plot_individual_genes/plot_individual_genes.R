here::i_am("plot_individual_genes/plot_individual_genes.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O ##
io$outdir <- file.path(io$basedir,"results/individual_genes")

## Define options ##

# Define cell types to plot
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
	"Blood_progenitors",
	# "Blood_progenitors_1",
	# "Blood_progenitors_2",
	"Erythroid",
	# "Erythroid1",
	# "Erythroid2",
	# "Erythroid3",
	"NMP",
	"Rostral_neurectoderm",
	"Caudal_neurectoderm",
	"Neural_crest",
	"Forebrain_Midbrain_Hindbrain",
	"Spinal_cord",
	"Surface_ectoderm"
	# "Visceral_endoderm",
	# "ExE_endoderm",
	# "ExE_ectoderm",
	# "Parietal_endoderm"
)

# Define classes to plot
# opts$classes <- c(
#   "E8.5_WT",
#   # "E8.5_Dnmt3aHET_Dnmt3bWT",
#   "E8.5_Dnmt3aKO_Dnmt3bWT",
#   "E8.5_Dnmt3aKO_Dnmt3bHET",
#   "E8.5_Dnmt3aWT_Dnmt3bKO",
#   "E8.5_Dnmt3aHET_Dnmt3bKO",
#   "E8.5_Dnmt3aKO_Dnmt3bKO",
#   "E8.5_Dnmt1KO"
# )
opts$classes <- c(
  # "E12.5_Dnmt3aWT_Dnmt3bHET",
  # "E12.5_Dnmt3b_KO",
  # "E12.5_Dnmt3a_HET_Dnmt3b_WT",
  # "E12.5_Dnmt3a_KO",
  # "Dnmt3a_HET_Dnmt3b_KO",
  # "Dnmt3a_HET_Dnmt3b_WT",
  # "Dnmt3a_KO_Dnmt3b_HET",
  "WT",
  "Dnmt3a_KO",
  "Dnmt3b_KO",
  "Dnmt1_KO",
  "Dnmt3ab_KO"
)

opts$rename_celltypes <- c(
  "Erythroid3" = "Erythroid",
  "Erythroid2" = "Erythroid",
  "Erythroid1" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors",
  "Anterior_Primitive_Streak" = "Primitive_Streak"
)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  .[,dataset:=ifelse(grepl("Grosswendt",sample),"CRISPR","KO")] %>% .[,dataset:=factor(dataset,levels=c("KO","CRISPR"))] %>%
  .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,opts$rename_celltypes)] %>% 
  .[pass_rnaQC==TRUE & celltype.mapped%in%opts$celltypes & class%in%opts$classes] %>%
  .[,class:=factor(class,levels=opts$classes)] %>% .[,celltype.mapped:=factor(celltype.mapped,levels=opts$celltypes)]

# Only consider cell types with sufficient observations in WT cells
celltypes.to.use <- sample_metadata[class=="WT",.(N=.N),by="celltype.mapped"] %>% .[N>=50,celltype.mapped]
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

################
## Parse data ##
################

# Remove underscores
# sample_metadata %>%
#   .[,class:=stringr::str_replace_all(class,"_"," ")] %>%
#   .[,class:=factor(class, levels=opts$classes %>% stringr::str_replace_all(.,"_"," "))] %>%
#   .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,"_"," ")] %>%
#   .[,celltype.mapped:=factor(celltype.mapped, levels=opts$celltypes %>% stringr::str_replace_all(.,"_"," "))]
  
# Rename celltypes
# .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,opts$rename_celltypes)] %>%
#   .[,celltype.mapped:=factor(celltype.mapped, levels=opts$celltypes)] 

##########
## Plot ##
##########

# genes.to.plot <- c("Tex19.1","Morc1","Dppa3","Rex1","Dppa5a","Dppa4","Dppa2","Zfp981")
# genes.to.plot <- c("Lefty1","Cd34","Tmsb4x","Fgf3","Spata7","Cer1","Spink1","Dppa4","Dppa5a","Prc1","Lefty2","Ube2c","Hba-x","Hbb-y","Hba-a1","Hbb-bh1")
# genes.to.plot <- c("Vegfa","Vegfb","Vegfc","Vegfd","Kdr","Flt1","Tal1","Runx1","Etv2)
# genes.to.plot <- c("Tet1","Tet2","Tet3","Dnmt1","Dnmt3a","Dnmt3b","Dnmt3l")
# genes.to.plot <- rownames(sce)[grep("Lefty",rownames(sce))]
# genes.to.plot <- fread(io$atlas.marker_genes) %>% .[celltype=="Epiblast" & score>=0.80,gene]
genes.to.plot <- fread(io$atlas.marker_genes)[,gene] %>% unique %>% head(n=3)
# genes.to.plot <- fread("/Users/ricard/data/gastrulation10x/results/differential/celltypes/E8.5/Neural_crest_vs_Forebrain_Midbrain_Hindbrain.txt.gz") %>% .[sig==T & logFC<(-2),gene]

for (i in 1:length(genes.to.plot)) {
  
  gene <- genes.to.plot[i]
  
  if (gene %in% rownames(sce)) {
    print(sprintf("%s/%s: %s",i,length(genes.to.plot),gene))
    outfile <- sprintf("%s/%s_boxplots_single_cells.pdf",io$outdir,gene)
  
    to.plot <- data.table(
      cell = colnames(sce),
      expr = logcounts(sce)[gene,]
    ) %>% merge(sample_metadata[,c("cell","sample","class","celltype.mapped","dataset")], by="cell") %>%
      .[,N:=.N,by=c("sample","celltype.mapped")] %>% .[N>=10]
    
    p_list <- list()
    for (j in levels(to.plot$dataset)) {
      p_list[[j]] <- ggplot(to.plot[dataset==j], aes(x=class, y=expr, fill=class)) +
        geom_violin(scale = "width", alpha=0.8) +
        geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.8) +
        # stat_summary(fun.data = give.n, geom = "text", size=3) +
        # scale_fill_brewer(palette="Dark2") +
        scale_fill_manual(values=opts$classes.colors) +
        facet_wrap(~celltype.mapped, scales="fixed") +
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
          legend.position = "none",
          legend.title = element_blank(),
          legend.text = element_text(size=rel(0.85))
        )
      }
      pdf(outfile, width=17, height=9)
      print(cowplot::plot_grid(plotlist=p_list, ncol = 2))
      dev.off()
  } else {
    print(sprintf("%s not found",gene))
  }
}

