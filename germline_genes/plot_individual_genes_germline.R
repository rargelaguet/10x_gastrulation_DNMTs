here::i_am("plot_individual_genes/plot_individual_genes.R")

source(here::here("settings.R"))
source(here::here("utils.R"))
source(here::here("plot_individual_genes/plot_utils.R"))

#####################
## Define settings ##
#####################

# I/O #
io$outdir <- file.path(io$basedir,"results/individual_genes/germline"); dir.create(io$outdir, showWarnings = F)

# Define classes to plot
opts$classes <- c(
  "WT", 
  "Dnmt3a_KO", 
  "Dnmt3b_KO",
  "Dnmt1_KO"
  # "Dnmt3ab_KO"
)

# Rename cell types
opts$rename_celltypes <- c(
  "Erythroid3" = "Erythroid",
  "Erythroid2" = "Erythroid",
  "Erythroid1" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors"
)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  .[pass_rnaQC==TRUE & celltype.mapped%in%opts$celltypes & class%in%opts$classes] %>%
  .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,opts$rename_celltypes)]

table(sample_metadata$class)
table(sample_metadata$celltype.mapped)

##########################
## Load germline genes ##
##########################

germline_genes <- c("Gtsf1", "Tex19.1", "Topaz1", "Rnf17", "Piwil2", "Mov10l1", "Asz1", "Ddx4", "Mael", "Fkbp6", "Tuba3b", "Sohlh2", "Rpl10l", "Dazl", "Asz1", "Hormad1")

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- load_SingleCellExperiment(io$sce, cells=sample_metadata$cell, normalise = TRUE)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

# Subset genes
genes.to.use <- intersect(germline_genes,rownames(sce))
sce <- sce[rownames(sce)%in%genes.to.use,]

##########
## Plot ##
##########

genes.to.plot <- genes.to.use
celltypes.to.plot <- c("Spinal_cord","Cardiomyocytes","Blood_progenitors","Gut","ExE_ectoderm", "ExE_endoderm")

for (i in genes.to.plot) {
  pdf(sprintf("%s/%s_germline.pdf",io$outdir,i), width=8, height=4)
  p <- plotting_fn(sce, gene = i, celltypes = celltypes.to.plot, classes = opts$classes) +
    facet_wrap(~celltype, scales="fixed", nrow=2)
  print(p)
  dev.off()
}
