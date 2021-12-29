here::i_am("plot_individual_genes/plot_individual_genes.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$outdir <- file.path(io$basedir,"results_new/individual_genes")

## options


##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  .[pass_rnaQC==TRUE & celltype.mapped%in%opts$celltypes & class%in%opts$classes]

table(sample_metadata$class)
table(sample_metadata$celltype.mapped)

# Subset cell types with enough observations
sample_metadata <- sample_metadata[,N:=.N,by="celltype.mapped"] %>% .[N>=25]

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- load_SingleCellExperiment(io$sce, cells=sample_metadata$cell, normalise = TRUE)
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

# Load Pluripotency genes
opts$pluripotency_genes <- fread("/Users/argelagr/data/gastrulation/rna/results/differential/E4.5Epiblast_vs_E7.5EctodermEndodermMesoderm.txt.gz") %>% 
  .[padj_fdr<=0.01 & logFC<(-3),symbol]
opts$pluripotency_genes <- opts$pluripotency_genes[opts$pluripotency_genes%in%rownames(sce)]

# Subset pluripotency genes
pluripotency.sce <- sce[opts$pluripotency_genes,]

#########
## PCA ##
#########



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
