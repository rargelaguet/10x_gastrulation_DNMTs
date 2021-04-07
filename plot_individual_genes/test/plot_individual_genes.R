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
	"Caudal_neurectoderm",
	"Neural_crest",
	"Forebrain_Midbrain_Hindbrain",
	"Spinal_cord",
	"Surface_ectoderm",
	"Visceral_endoderm",
	"ExE_endoderm",
	"ExE_ectoderm",
	"Parietal_endoderm"
)

# Define classes to plot
opts$classes <- c(
  "E8.5_WT",
  "E8.5_Dnmt3aKO_Dnmt3bWT",
  "E8.5_Dnmt3aHET_Dnmt3bKO",
  # "E8.5_Dnmt3aHET_Dnmt3bWT",
  # "E8.5_Dnmt3aKO_Dnmt3bHET",
  "E8.5_Dnmt3aWT_Dnmt3bKO",
  "E8.5_Dnmt3aKO_Dnmt3bKO",
  "E8.5_Dnmt1KO"
)

# Update sample metadata
sample_metadata <- fread(io$metadata) %>% 
  .[pass_QC==TRUE & class%in%opts$classes & celltype.mapped%in%opts$celltypes] %>%
  .[,celltype.mapped:=factor(celltype.mapped, levels=opts$celltypes)] %>%
  .[,class:=factor(class,levels=opts$classes)]

opts$min.cells <- 25
sample_metadata <- sample_metadata %>%
  .[,N:=.N,by=c("celltype.mapped","sample")] %>% .[N>=opts$min.cells] %>% .[,N:=NULL] %>%
  droplevels


table(sample_metadata$class)
table(sample_metadata$sample)
table(sample_metadata$celltype.mapped)

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- load_SingleCellExperiment(io$sce, cells=sample_metadata$cell, normalise = FALSE)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

# Normalise
sce <- logNormCounts(sce)

##########
## Plot ##
##########

# genes.to.plot <- fread("/Users/ricard/data/gastrulation10x/results/differential/celltypes/E8.5/Neural_crest_vs_Forebrain_Midbrain_Hindbrain.txt.gz") %>%
  # .[sig==T & logFC<0,gene]

# genes.to.plot <- fread("/Users/ricard/data/10x_gastrulation_DNMTs/results/differential/E8.5_Dnmt3aKO_Dnmt3bKO_vs_E8.5_WT_Epiblast.txt.gz") %>%
#   setorder(-log_padj_fdr,na.last=T) %>%
  # .[sig==T & logFC>0 & abs(logFC)>2,gene] %>% head(n=15)

genes.to.plot <- fread(io$atlas.marker_genes)$gene %>% unique

# genes.to.plot <- rownames(sce)[grep("^Xlr",rownames(sce))]
# genes.to.plot <- rownames(sce)[grep("^Rhox",rownames(sce))]
# genes.to.plot <- rownames(sce)[grep("^Tet",rownames(sce))]
# genes.to.plot <- rownames(sce)[grep("^Dnmt",rownames(sce))]

# genes.to.plot <- genes.to.plot[genes.to.plot%in%rownames(sce)]

for (i in 1:length(genes.to.plot)) {
  
  gene <- genes.to.plot[i]
  
  if (gene %in% rownames(sce)) {
    print(sprintf("%s/%s: %s",i,length(genes.to.plot),gene))
    outfile <- sprintf("%s/%s.pdf",io$outdir,gene)
    
    if (!file.exists(outfile)) {
      
      to.plot <- data.table(
        cell = colnames(sce),
        expr = logcounts(sce)[gene,]
      ) %>% merge(sample_metadata[,c("cell","class","celltype.mapped")], by="cell")
      
      p <- ggplot(to.plot, aes(x=class, y=expr, fill=class)) +
        geom_violin(scale = "width", alpha=0.8) +
        geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.8) +
        stat_summary(fun.data = give.n, geom = "text", size=2.5) + 
        facet_wrap(~celltype.mapped, scales="fixed") +
        theme_classic() +
        labs(title=gene, x="",y=sprintf("%s expression",gene)) +
        guides(x = guide_axis(angle = 90)) +
        theme(
          strip.text = element_text(size=rel(0.7)),
          # plot.title = element_text(hjust = 0.5, size=rel(1.1), color="black"),
          plot.title = element_blank(),
          axis.text.x = element_text(colour="black",size=rel(0.7)),
          # axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(colour="black",size=rel(1.0)),
          axis.title.y = element_text(colour="black",size=rel(1.0)),
          legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(size=rel(0.75))
        )
          
        pdf(outfile, width=11, height=10)
        # jpeg(sprintf("%s/%s.jpeg",io$outdir,gene), width = 1000, height = 800)
        print(p)
        dev.off()
        
    } else {
      print(sprintf("%s already exists...",outfile))
    }

  } else {
    print(sprintf("%s not found",gene))
  }
}
