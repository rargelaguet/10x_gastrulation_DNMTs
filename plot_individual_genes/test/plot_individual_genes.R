suppressPackageStartupMessages(library(SingleCellExperiment))

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
  source("/Users/ricard/10x_gastrulation_DNMTs/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/10x_gastrulation_DNMTs/settings.R")
  source("/homes/ricard/10x_gastrulation_DNMTs/utils.R")
}

io$outdir <- paste0(io$basedir,"/results/individual_genes/test")

# Define cell types to plot
opts$celltypes <- c(
	# "Epiblast",
	# "Primitive_Streak",
	# "Caudal_epiblast",
	# "PGC",
	# "Anterior_Primitive_Streak",
	# "Notochord",
	# "Def._endoderm",
	# "Gut",
	# "Nascent_mesoderm",
	# "Mixed_mesoderm",
	# "Intermediate_mesoderm",
	# "Caudal_Mesoderm",
	# "Paraxial_mesoderm",
	# "Somitic_mesoderm",
	# "Pharyngeal_mesoderm",
	# "Cardiomyocytes",
	# "Allantois",
	# "ExE_mesoderm",
	# "Mesenchyme",
	# "Haematoendothelial_progenitors",
	# "Endothelium",
	# "Blood_progenitors_1",
	# "Blood_progenitors_2",
	# "Erythroid1",
	# "Erythroid2",
	# "Erythroid3",
	# "NMP",
	# "Rostral_neurectoderm",
	# "Caudal_neurectoderm",
	"Neural_crest"
	# "Forebrain_Midbrain_Hindbrain",
	# "Spinal_cord",
	# "Surface_ectoderm",
	# "Visceral_endoderm",
	# "ExE_endoderm",
	# "ExE_ectoderm"
	# "Parietal_endoderm"
)

# Define classes to plot
opts$classes <- c(
  "E8.5_WT",
  "E8.5_Dnmt3aKO_Dnmt3bWT",
  "E8.5_Dnmt3aHET_Dnmt3bKO",
  # "E8.5_Dnmt3aHET_Dnmt3bWT",
  # "E8.5_Dnmt3aKO_Dnmt3bHET",
  "E8.5_Dnmt3aWT_Dnmt3bKO",
  "E8.5_Dnmt3aKO_Dnmt3bKO"
)

# Update sample metadata
sample_metadata <- fread(io$metadata) %>% 
  .[pass_QC==TRUE & class%in%opts$classes & celltype.mapped%in%opts$celltypes] %>%
  # .[batch%in%opts$batches & celltype.mapped%in%opts$celltypes] %>%
  .[,celltype.mapped:=factor(celltype.mapped, levels=opts$celltypes)] %>%
  .[,class:=factor(class,levels=opts$classes)]

# Rename samples
foo <- sample_metadata[,c("batch","class")] %>% unique %>% .[,sample:=paste(class,1:.N,sep="_"),by="class"]
sample_metadata <- sample_metadata %>% merge(foo,by=c("batch","class"))

table(sample_metadata$class)
table(sample_metadata$sample)
table(sample_metadata$celltype.mapped)

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- load_SingleCellExperiment(io$sce, cells=sample_metadata$cell, normalise = TRUE, remove_non_expressed_genes = TRUE)

##########
## Plot ##
##########

genes.to.plot <- fread("/Users/ricard/data/gastrulation10x/results/differential/celltypes/E8.5/Neural_crest_vs_Forebrain_Midbrain_Hindbrain.txt.gz") %>%
  .[sig==T & logFC<0,gene]

# genes.to.plot <- rownames(sce)[grep("Tet",rownames(sce))]
# genes.to.plot <- rownames(sce)[grep("Dnmt",rownames(sce))]

genes.to.plot <- genes.to.plot[genes.to.plot%in%rownames(sce)]

for (i in 1:length(genes.to.plot)) {
  gene <- genes.to.plot[i]
  print(sprintf("%s/%s: %s",i,length(genes.to.plot),gene))
  
  # Create data.table to plot
  to.plot <- data.table(
    cell = colnames(sce),
    expr = logcounts(sce[gene,])[1,]
  ) %>% merge(sample_metadata, by="cell")


  p <- ggplot(to.plot, aes(x=class, y=expr)) +
    geom_violin(scale = "width", alpha=0.8, fill="grey70") +
    geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.8, fill="grey70") +
    stat_summary(fun.data = give.n, geom = "text") + 
    # facet_wrap(~celltype.mapped, scales="fixed") +
    theme_classic() +
    labs(title=gene, x="",y=sprintf("%s expression",gene)) +
    theme(
      # strip.text = element_text(size=rel(1.0)),
      # plot.title = element_text(hjust = 0.5, size=rel(1.1), color="black"),
      plot.title = element_blank(),
      axis.text.x = element_text(colour="black",size=rel(1.1)),
      # axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(colour="black",size=rel(1.0)),
      axis.title.y = element_text(colour="black",size=rel(1.2)),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size=rel(1.1))
    )
    
  # pdf(sprintf("%s/%s.pdf",io$outdir,i), width=8, height=3.5, useDingbats = F)
  # ggsave("ggtest.png", width = 3.25, height = 3.25, dpi = 1200)
  jpeg(sprintf("%s/%s.jpeg",io$outdir,gene), width = 800, height = 400)
  print(p)
  dev.off()
}
