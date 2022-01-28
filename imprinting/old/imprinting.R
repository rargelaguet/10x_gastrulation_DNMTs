source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

io$imprinted.genes <- "/Users/argelagr/data/mm10_regulation/imprinting/parsed/mousebook_imprinted_genes.txt.gz"
io$outdir <- paste0(io$basedir,"/results_all/imprinting"); dir.create(io$outdir, showWarnings = F)

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
	# "Blood_progenitors",
	"Blood_progenitors_1",
	"Blood_progenitors_2",
	# "Erythroid",
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
	"ExE_ectoderm"
	# "Parietal_endoderm"
)

# Define classes to plot
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

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  .[pass_rnaQC==TRUE & class%in%opts$classes & celltype.mapped%in%opts$celltypes] %>%
  .[,celltype.mapped:=factor(celltype.mapped, levels=opts$celltypes)] %>%
  .[,class:=factor(class,levels=opts$classes)]

table(sample_metadata$class)
table(sample_metadata$sample)
table(sample_metadata$celltype.mapped)

##########################
## Load imprinted genes ##
##########################

imprinting.dt <- fread(io$imprinted.genes) %>%
  setnames(c("gene","allele"))

##############################
## Load RNA expression data ##
##############################

# Load SingleCellExperiment object
sce <- load_SingleCellExperiment(
  file = io$sce, 
  features = imprinting.dt$gene,
  cells = sample_metadata$cell, 
  normalise = TRUE
)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

################
## Parse data ##
################

sce <- sce[rowSums(counts(sce))>100,]

to.plot <- data.table(
  cell = colnames(sce),
  imprinting_expr = colMeans(logcounts(sce))
) %>% merge(sample_metadata,by="cell")

##########
## Plot ##
##########

p <- ggplot(to.plot, aes(x=class, y=imprinting_expr, fill=class)) +
  geom_violin(scale = "width", alpha=0.8) +
  geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.8) +
  scale_x_discrete(drop=F) +
  theme_classic() +
  labs(x="",y="Imprinting expression") +
  guides(x = guide_axis(angle = 90)) +
  theme(
    plot.title = element_blank(),
    axis.text.x = element_text(colour="black",size=rel(0.95)),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(colour="black",size=rel(1.0)),
    axis.title.y = element_text(colour="black",size=rel(1.0)),
    legend.position = "none"
  )


pdf(paste0(io$outdir,"/imprinting_expr.pdf"), width=7, height=4)
print(p)
dev.off()

p <- ggplot(to.plot, aes(x=class, y=imprinting_expr, fill=class)) +
  geom_violin(scale = "width", alpha=0.8) +
  geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.8) +
  # geom_jitter(size=2, shape=21, stroke=0.2, alpha=0.5) +
  # scale_fill_manual(values=opts$classes.colors, drop=F) +
  scale_x_discrete(drop=F) +
  stat_summary(fun.data = give.n, geom = "text", size=2.5) +
  facet_wrap(~celltype.mapped, scales="fixed") +
  theme_classic() +
  labs(x="",y="Imprinting expression") +
  guides(x = guide_axis(angle = 90)) +
  theme(
    strip.text = element_text(size=rel(0.85)),
    # plot.title = element_text(hjust = 0.5, size=rel(1.1), color="black"),
    plot.title = element_blank(),
    axis.text.x = element_text(colour="black",size=rel(0.95)),
    # axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(colour="black",size=rel(1.0)),
    axis.title.y = element_text(colour="black",size=rel(1.0)),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size=rel(0.85))
  )


pdf(paste0(io$outdir,"/imprinting_expr_by_celltype.pdf"), width=11, height=10)
print(p)
dev.off()