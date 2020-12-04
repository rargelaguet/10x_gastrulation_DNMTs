suppressPackageStartupMessages(library(SingleCellExperiment))

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/10x_gastrulation_DNMTs/settings.R")
}

io$outdir <- paste0(io$basedir,"/results/individual_genes")

# Define cell types to plot
opts$celltypes <- c(
	# "Epiblast",
	# "Primitive_Streak",
	# "Caudal_epiblast",
	# "PGC",
	# "Anterior_Primitive_Streak",
	"Notochord",
	"Def._endoderm",
	"Gut",
	# "Nascent_mesoderm",
	# "Mixed_mesoderm",
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
	# "Visceral_endoderm",
	"ExE_endoderm",
	"ExE_ectoderm"
	# "Parietal_endoderm"
)

# Define classes to plot
opts$classes <- c(
  # "E12.5_Dnmt3aWT_Dnmt3bHET",
  # "E12.5_Dnmt3aWT_Dnmt3bKO",
  # "E12.5_Dnmt3aHET_Dnmt3bWT",
  # "E12.5_Dnmt3aKO_Dnmt3bWT",
  "E8.5_Dnmt3aKO_Dnmt3bWT",
  "E8.5_Dnmt3aWT_Dnmt3bWT",
  # "E8.5_Dnmt3aHET_Dnmt3bKO",
  # "E8.5_Dnmt3aHET_Dnmt3bWT",
  # "E8.5_Dnmt3aKO_Dnmt3bHET",
  "E8.5_Dnmt3aKO_Dnmt3bKO",
  "E8.5_Dnmt3aWT_Dnmt3bKO"
)

# Define batches to plot
# opts$batches <- c(
#   "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001", 
#   "SIGAB6_E85_3_Dnmt3aWT_Dnmt3b_WT_L002", 
#   # "SIGAC6_E85_5_Dnmt3aKO_Dnmt3b_Het_L003", 
#   # "SIGAD6_E85_8_Dnmt3aHet_Dnmt3b_KO_L004",
#   "15_E8_5_D3A_WT_D3B_WT_L007",
#   "17_E8_5_D3A_KO_D3B_WT_L008",
#   "2_E8_5_D3A_WT_D3B_KO_L003",
#   # "3_E8_5_D3A_HET_D3B_WT_L004",
#   "7_E8_5_D3A_WT_D3B_KO_L005",
#   "8_E8_5_D3A_KO_D3B_KO_L006"
# )

# Update sample metadata
sample_metadata <- sample_metadata %>% 
  .[class%in%opts$classes & celltype.mapped%in%opts$celltypes] %>%
  # .[batch%in%opts$batches & celltype.mapped%in%opts$celltypes] %>%
  .[,celltype.mapped:=factor(celltype.mapped, levels=opts$celltypes)]

table(sample_metadata$class)
table(sample_metadata$batch)
table(sample_metadata$celltype.mapped)

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- readRDS(io$sce)[,sample_metadata$cell]

# Remove genes that are not expressed
sce <- sce[rowMeans(counts(sce))>0,]

# Load gene metadata
# gene_metadata <- fread(io$gene_metadata) %>%
#   .[symbol!="" & ens_id%in%rownames(sce)]

################
## Parse data ##
################

## Normalisation
# sce <- batchelor::multiBatchNorm(sce, batch=as.factor(sce$batch))
# sce <- scater::logNormCounts(sce)

# Rename genes
# new.names <- gene_metadata$symbol
# names(new.names) <- gene_metadata$ens_id
# sce <- sce[rownames(sce) %in% names(new.names),]
# rownames(sce) <- new.names[rownames(sce)]
# stopifnot(sum(is.na(new.names))==0)
# stopifnot(sum(duplicated(new.names))==0)

##########
## Plot ##
##########

# genes.to.plot <- c("Cd59a", "Cdx2", "Gata4", "Cd9")
# genes.to.plot <- rownames(sce)
# genes.to.plot <- c("Dppa4","Erdr1","Slc25a31","Uba52","Pim2","Slc7a3", "Tex19.1","Cdx1","Fosb","Jun","Fos","Pou5f1","Pim2","Slc7a3") %>% unique
# genes.to.plot <- rownames(sce)[grep("Psm",rownames(sce))]
# genes.to.plot <- rownames(sce)[grep("Tet",rownames(sce))]
genes.to.plot <- rownames(sce)[grep("Dnmt",rownames(sce))]
# genes.to.plot <- fread("/Users/ricard/data/mm10_regulation/imprinting/mousebook_imprinted_genes.txt.gz") %>% .[symbol%in%rownames(sce),symbol] 

for (i in 1:length(genes.to.plot)) {
  gene <- genes.to.plot[i]
  print(sprintf("%s/%s: %s",i,length(genes.to.plot),gene))
  
  # Create data.table to plot
  to.plot <- data.table(
    cell = colnames(sce),
    expr = logcounts(sce[gene,])[1,]
  ) %>% merge(sample_metadata, by="cell")

  # Plot
  # p <- ggplot(to.plot, aes(x=celltype.mapped, y=expr, fill=celltype.mapped)) +
  #   geom_violin(scale = "width") +
  #   geom_boxplot(width=0.5, outlier.shape=NA) +
  #   scale_fill_manual(values=opts$celltype.colors) +
  #   facet_wrap(~class) +
  #   theme_classic() +
  #   labs(x="",y="RNA expression") +
  #   theme(
  #     axis.text.x = element_text(colour="black",size=rel(1.0), angle=50, hjust=1),
  #     axis.text.y = element_text(colour="black",size=rel(1.0)),
  #     legend.position="none"
  #   )
  
  p <- ggplot(to.plot, aes(x=class, y=expr, fill=class)) +
    geom_violin(scale = "width", alpha=0.8) +
    geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.8) +
    # geom_jitter(size=2, shape=21, stroke=0.2, alpha=0.5) +
    # scale_fill_manual(values=opts$colors) +
    facet_wrap(~celltype.mapped, scales="fixed") +
    theme_classic() +
    labs(title=gene, x="",y=sprintf("%s expression",gene)) +
    theme(
      strip.text = element_text(size=rel(1.0)),
      # plot.title = element_text(hjust = 0.5, size=rel(1.1), color="black"),
      plot.title = element_blank(),
      # axis.text.x = element_text(colour="black",size=rel(1.2), angle=50, hjust=1),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(colour="black",size=rel(1.0)),
      axis.title.y = element_text(colour="black",size=rel(1.2)),
      legend.position="top",
      legend.title = element_blank(),
      legend.text = element_text(size=rel(1.1))
    )
    
  # pdf(sprintf("%s/%s.pdf",io$outdir,i), width=8, height=3.5, useDingbats = F)
  # ggsave("ggtest.png", width = 3.25, height = 3.25, dpi = 1200)
  jpeg(sprintf("%s/%s.jpeg",io$outdir,gene), width = 1300, height = 800)
  print(p)
  dev.off()
}
