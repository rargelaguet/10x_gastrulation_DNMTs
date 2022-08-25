here::i_am("plot_individual_genes/plot_individual_genes.R")

source(here::here("settings.R"))
source(here::here("utils.R"))
source(here::here("plot_individual_genes/plot_utils.R"))
source(here::here("repeats/utils.R"))

#####################
## Define settings ##
#####################

# I/O #
io$outdir <- file.path(io$basedir,"results/individual_genes/repeats"); dir.create(io$outdir, showWarnings = F, recursive = T)

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


###############
## Load data ##
###############

exp <- fread(paste0(io$basedir, "/repeats/repeats_expr.txt.gz")) %>% 
  .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,opts$rename_celltypes)]



# Subset genes
# genes.to.use <- c(
#   "LINE_L1",
#   "IAP",
#   "major_satellite",
#   "minor_satellite",
#   "LTR_ERV1"
# )

genes.to.use <- exp[, unique(repeat_class)]


##########
## Plot ##
##########

genes.to.plot <- genes.to.use
# celltypes.to.plot <- c(
#   "Spinal_cord",
#   "Cardiomyocytes",
#   "Blood_progenitors",
#   "Gut",
#   "ExE_ectoderm", 
#   "ExE_endoderm"
# )

# which celltypes have good coverage in Dnmt1 KO?
sample_metadata[, .N, .(celltype.mapped, dataset, class)][N>50][class=="Dnmt1_KO" & dataset=="KO"][, celltype.mapped]

celltypes.to.plot <- c(
  "ExE_ectoderm",
  "Erythroid",
  "Rostral_neurectoderm",
  "ExE_endoderm",
  "Gut",
  "Surface_ectoderm"
)

for (i in genes.to.plot) {
  pdf(sprintf("%s/%s_repeats.pdf",io$outdir,i), width=8, height=4)
  p <- plotting_fn2(exp, gene = i, celltypes = celltypes.to.plot, classes = opts$classes) +
    facet_wrap(~celltype, scales="fixed", nrow=2)
  print(p)
  dev.off()
}
