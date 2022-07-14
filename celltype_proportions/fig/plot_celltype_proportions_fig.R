here::i_am("celltype_proportions/plot_celltype_proportions.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

io$metadata <- file.path(io$basedir,"results/mapping/sample_metadata_after_mapping.txt.gz")
io$outdir <- file.path(io$basedir,"results/celltype_proportions/fig")

# I/O
dir.create(io$outdir, showWarnings = F)

# Options
opts$remove_ExE_cells <- FALSE

opts$celltypes = c(
  "Epiblast",
  # "Primitive_Streak",
  "Caudal_epiblast",
  # "PGC",
  # "Anterior_Primitive_Streak",
  # "Notochord",
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
  # "Haematoendothelial_progenitors",
  "Haematoend._progenitors",
  "Endothelium",
  "Blood_progenitors",
  # "Blood_progenitors_1",
  # "Blood_progenitors_2",
  # "Erythroid1",
  # "Erythroid2",
  # "Erythroid3",
  "Erythroid",
  "NMP",
  "Rostral_neurectoderm",
  "Caudal_neurectoderm",
  "Neural_crest",
  "Brain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  # "ExE_endoderm",
  "ExE_ectoderm"
  # "Parietal_endoderm"
)

opts$rename_celltypes <- c(
  "Erythroid3" = "Erythroid",
  "Erythroid2" = "Erythroid",
  "Erythroid1" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors",
  "Forebrain_Midbrain_Hindbrain" = "Brain",
  "Haematoendothelial_progenitors" = "Haematoend._progenitors"
  # "Intermediate_mesoderm" = "Mixed_mesoderm",
  # "Paraxial_mesoderm" = "Mixed_mesoderm",
  # "Nascent_mesoderm" = "Mixed_mesoderm",
  # "Pharyngeal_mesoderm" = "Mixed_mesoderm"
  # "Visceral_endoderm" = "ExE_endoderm"
)

opts$classes <- c("WT", "Dnmt3a_KO", "Dnmt3b_KO", "Dnmt1_KO")

# opts$samples <- c(
#   "WT_1", "WT_2", "WT_3", "WT_4", "WT_6",
#   "Dnmt3a_KO_1", "Dnmt3a_KO_12", "Dnmt3a_KO_2", "Dnmt3a_KO_13", "Dnmt3a_KO_14",
#   "Dnmt3b_KO_1", "Dnmt3b_KO_2", "Dnmt3b_KO_6", "Dnmt3b_KO_7", "Dnmt3b_KO_9", 
#   "Dnmt1_KO_10", "Dnmt1_KO_15", "Dnmt1_KO_2", "Dnmt1_KO_3", "Dnmt1_KO_9"
#   # "Dnmt3ab_KO_1", "Dnmt3ab_KO_2"
# )

##########################
## Load sample metadata ##
##########################

cell_metadata.dt <- fread(io$metadata) %>%
  setnames("celltype.mapped","celltype") %>%
  .[,celltype:=stringr::str_replace_all(celltype,opts$rename_celltypes)] %>%
  .[pass_rnaQC==TRUE & class%in%opts$classes & celltype%in%opts$celltypes] %>%
  .[,dataset:=ifelse(grepl("Grosswendt",sample),"CRISPR","KO")] %>%
  .[,c("cell","sample","alias","class","celltype","dataset")]

if (opts$remove_ExE_cells) {
  print("Removing ExE cells...")
  cell_metadata.dt <- cell_metadata.dt %>%
    .[!celltype%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

###############
## Per class ##
###############

to.plot <- cell_metadata.dt %>%
  .[,N:=.N,by="class"] %>%
  .[,.(N=.N, celltype_proportion=.N/unique(N)),by=c("class","alias","celltype")] 

# Define colours and cell type order
opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors) %in% unique(to.plot$celltype)]
to.plot[,celltype:=factor(celltype, levels=rev(names(opts$celltype.colors)))]

p <- ggplot(to.plot, aes(x=celltype, y=N)) +
  geom_bar(aes(fill=celltype), stat="identity", color="black",) +
  scale_fill_manual(values=opts$celltype.colors, drop=F) +
  facet_wrap(~alias, nrow=4, scales="free_x") +
  coord_flip() +
  labs(y="Number of cells") +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(color="black", size=rel(0.8)),
    axis.title.x = element_text(color="black", size=rel(0.9)),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size=rel(0.75), color="black"),
    axis.text.x = element_text(size=rel(1), color="black")
  )
  
pdf(file.path(io$outdir,"celltype_proportions_horizontal_barplots.pdf"), width=10, height=14)
print(p)
dev.off()


############################
## Per class and data set ##
############################

# i <- "WT"
for (i in opts$classes) {
  
  to.plot <- cell_metadata.dt %>% 
    .[class==i] %>%
    .[,N:=.N,by=c("class","alias")] %>%
    .[,.(N=.N, celltype_proportion=.N/unique(N)),by=c("class","alias","celltype","dataset")] 
  
  # unique(to.plot$celltype)[!unique(to.plot$celltype)%in%opts$celltypes]
  # celltype.order <- opts$celltypes[opts$celltypes%in%unique(to.plot$celltype)] %>% rev
  celltype.order <- to.plot %>% .[,mean(celltype_proportion),by="celltype"] %>% setorder(-V1) %>% .$celltype
  to.plot <- to.plot %>% .[,celltype:=factor(celltype,levels=celltype.order)]
  
  p <- ggplot(to.plot, aes(x=celltype, y=celltype_proportion, fill=celltype)) +
    geom_point(shape=21, size=1.25, stroke=0.1) +
    geom_boxplot(alpha=0.75, outlier.shape=NA) +
    facet_wrap(~dataset, scales = "free_x", nrow=1) +
    coord_flip() +
    scale_fill_manual(values=opts$celltype.colors[names(opts$celltype.colors)%in%unique(to.plot$celltype)], drop=F) +
    theme_classic() +
    labs(y="Number of cells", x="", title=i) +
    theme(
      legend.position = "none",
      # axis.title = element_blank(),
      plot.title = element_text(size=rel(1.25), hjust=0.5, color="black"),
      axis.text.y = element_text(color="black"),
      axis.text.x = element_text(color="black")
    )
  
  pdf(file.path(io$outdir,sprintf("%s_boxplots_celltype_proportions_per_dataset.pdf",i)), width=9, height=5)
  print(p)
  dev.off()
}
