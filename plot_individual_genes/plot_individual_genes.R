suppressPackageStartupMessages(library(SingleCellExperiment))

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/10x_gastrulation_DNMTs/settings.R")
}

## I/O ##

io$outdir <- paste0(io$basedir,"/results/individual_genes")

## Define options ##

# Define cell types to plot
opts$celltypes = c(
	"Epiblast",
	"Primitive_Streak",
	"Caudal_epiblast",
	"PGC",
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
  # "E8.5_Dnmt3aKO_Dnmt3bWT", 
  "E8.5_Dnmt3aWT_Dnmt3bWT", 
  # "E8.5_Dnmt3aKO_Dnmt3bHET", 
  "E8.5_Dnmt3aHET_Dnmt3bKO"
  # "E12.5_Dnmt3aWT_Dnmt3bHET", 
  # "E12.5_Dnmt3aWT_Dnmt3bKO" 
)

# Update sample metadata
sample_metadata <- sample_metadata %>% 
  .[class%in%opts$classes & celltype.mapped%in%opts$celltypes] %>%
  .[,celltype.mapped:=factor(celltype.mapped, levels=opts$celltypes)]

table(sample_metadata$class)
table(sample_metadata$celltype.mapped)

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- readRDS(io$sce)[,sample_metadata$cell]

# Remove genes that are not expressed
sce <- sce[rowMeans(counts(sce))>0,]

# Load gene metadata
gene_metadata <- fread(io$gene_metadata) %>%
  .[symbol!="" & ens_id%in%rownames(sce)]

################
## Parse data ##
################

# Rename genes
new.names <- gene_metadata$symbol
names(new.names) <- gene_metadata$ens_id
sce <- sce[rownames(sce) %in% names(new.names),]
rownames(sce) <- new.names[rownames(sce)]

# Sanity cehcks
stopifnot(sum(is.na(new.names))==0)
stopifnot(sum(duplicated(new.names))==0)

##########
## Plot ##
##########

# genes.to.plot <- c("Dnmt3l","Apoe")
# genes.to.plot <- rownames(sce)
genes.to.plot <- "Xist"

for (i in 1:length(genes.to.plot)) {
  gene <- genes.to.plot[i]
  print(sprintf("%s/%s: %s",i,length(genes.to.plot),gene))
  
  # Create data.table to plot
  to.plot <- data.table(
    cell = colnames(sce),
    expr = logcounts(sce[gene,])[1,]
  ) %>% merge(sample_metadata, by="cell")

  # Plot
  p <- ggplot(to.plot, aes(x=celltype.mapped, y=expr, fill=celltype.mapped)) +
    # geom_jitter(size=0.8) +
    geom_violin(scale = "width") +
    geom_boxplot(width=0.5, outlier.shape=NA) +
    scale_fill_manual(values=opts$celltype.colors) +
    facet_wrap(~class) +
    theme_classic() +
    labs(x="",y="RNA expression") +
    theme(
      axis.text.x = element_text(colour="black",size=rel(1.2), angle=50, hjust=1),
      axis.text.y = element_text(colour="black",size=rel(1.0)),
      legend.position="none"
    )
    
  # pdf(sprintf("%s/%s.pdf",io$outdir,i), width=8, height=3.5, useDingbats = F)
  # ggsave("ggtest.png", width = 3.25, height = 3.25, dpi = 1200)
  jpeg(sprintf("%s/%s.jpeg",io$outdir,gene), width = 1300, height = 400)
  print(p)
  dev.off()
}
