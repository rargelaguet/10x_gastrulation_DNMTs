here::i_am("pluripotency_score/pluripotency_score.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$outdir <- file.path(io$basedir,"results_new/pluripotency_score"); dir.create(io$outdir)

# Options
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
  # "Mixed_mesoderm",
  "Intermediate_mesoderm",
  # "Caudal_Mesoderm",
  "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "Pharyngeal_mesoderm",
  "Cardiomyocytes",
  "Allantois",
  # "ExE_mesoderm",
  # "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  # "Blood_progenitors",
  # "Blood_progenitors_1",
  "Blood_progenitors_2",
  # "Erythroid",
  # "Erythroid1",
  # "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord"
  # "Surface_ectoderm"
  # "Visceral_endoderm",
  # "ExE_endoderm",
  # "ExE_ectoderm",
  # "Parietal_endoderm"
)

opts$classes <- c("E8.5_WT", "E8.5_Dnmt3aKO_Dnmt3bWT", "E8.5_Dnmt3aWT_Dnmt3bKO", "E8.5_Dnmt3aKO_Dnmt3bKO", "E8.5_Dnmt1KO")

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  .[pass_rnaQC==TRUE & celltype.mapped%in%opts$celltypes]
  # .[pass_rnaQC==TRUE & celltype.mapped%in%opts$celltypes & class%in%opts$classes]

table(sample_metadata$class)
table(sample_metadata$celltype.mapped)

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- load_SingleCellExperiment(io$sce, cells=sample_metadata$cell, normalise = TRUE)
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

# Load Pluripotency genes


# Subset pluripotency genes


#######################
## Load marker genes ##
#######################

# opts$pluripotency_genes <- fread("/Users/argelagr/data/gastrulation/rna/results/differential/E4.5Epiblast_vs_E7.5EctodermEndodermMesoderm.txt.gz") %>% 
#   .[padj_fdr<=0.01 & logFC<(-3),symbol]
# opts$pluripotency_genes <- opts$pluripotency_genes[opts$pluripotency_genes%in%rownames(sce)]

pluripotency.genes <- fread(io$atlas.marker_genes) %>% .[celltype=="Epiblast" & score>=0.80,gene]

#################
## Filter data ##
#################

pluripotency.sce <- sce[pluripotency.genes,]

tmp <- pluripotency.sce[,pluripotency.sce$class%in%c("E8.5_WT","E8.5_Dnmt1KO","E8.5_Dnmt3aKO_Dnmt3bKO")]

#########
## PCA ##
#########

pca <- prcomp(t(logcounts(tmp)))

# pca_features.dt <- pca$rotation %>%
#   as.data.table(keep.rownames = T) %>%
#   setnames("rn","gene")

pca.dt <- pca$x %>% 
  as.data.table(keep.rownames = T) %>%
  setnames("rn","cell")

celltypes.to.plot <- c("Caudal_epiblast","Caudal_neurectoderm","Rostral_neurectoderm","Spinal_cord","Gut","Cardiomyocytes")

to.plot <- pca.dt %>%
  merge(sample_metadata[celltype.mapped%in%celltypes.to.plot],by="cell") %>% .[,celltype.mapped:=factor(celltype.mapped,levels=celltypes.to.plot)] %>%
  .[,N:=.N,by=c("celltype.mapped","class")] %>% .[N>=30] %>% .[,N:=NULL] %>%
  # .[,N:=length(unique(class)),by="celltype.mapped"] %>% .[N>=3] %>% .[,N:=NULL] %>%
  .[,pluripotency_score:=abs(PC1)]

p <- ggplot(to.plot, aes(x=class, y=pluripotency_score, fill=class)) +
  # geom_jitter(size=0.8) +
  geom_violin(scale = "width") +
  geom_boxplot(width=0.5, outlier.shape=NA) +
  facet_wrap(~celltype.mapped) +
  # geom_hline(yintercept=0) +
  scale_fill_brewer(palette="Dark2") +
  # scale_fill_manual(values=opts$celltype.colors) +
  theme_classic() +
  labs(x="",y="Pluripotency score (PC1 of Epiblast marker genes)") +
  theme(
    # axis.text.x = element_text(colour="black",size=rel(1.0), angle=50, hjust=1),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(colour="black",size=rel(1.0)),
    legend.position="top",
    legend.title=element_blank()
  )

# pdf(sprintf("%s/%s.pdf",io$outdir,i), width=8, height=3.5, useDingbats = F)
print(p)
# dev.off()

##########
## Plot ##
##########

to.plot <- data.table(
  cell = colnames(pluripotency.sce),
  pluripotency_score = colMeans(logcounts(pluripotency.sce))
) %>% merge(sample_metadata,by="cell")

p <- ggplot(to.plot, aes(x=class, y=pluripotency_score, fill=alias)) +
  # geom_jitter(size=0.8) +
  geom_violin(scale = "width") +
  geom_boxplot(width=0.5, outlier.shape=NA) +
  # scale_fill_manual(values=opts$celltype.colors) +
  facet_wrap(~celltype.mapped) +
  theme_classic() +
  labs(x="",y="RNA expression") +
  theme(
    # axis.text.x = element_text(colour="black",size=rel(1.0), angle=50, hjust=1),
    axis.text.x = element_blank(),
    axis.text.y = element_text(colour="black",size=rel(1.0)),
    legend.position="none",
    legend.title=element_blank()
  )

# pdf(sprintf("%s/%s.pdf",io$outdir,i), width=8, height=3.5, useDingbats = F)
# ggsave("ggtest.png", width = 3.25, height = 3.25, dpi = 1200)
# jpeg(sprintf("%s/%s.jpeg",io$outdir,gene), width = 1500, height = 600)
print(p)
# dev.off()
