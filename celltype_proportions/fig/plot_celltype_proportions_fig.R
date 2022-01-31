here::i_am("celltype_proportions/plot_celltype_proportions.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
# p$add_argument('--samples',         type="character",       nargs="+",   help='Samples')
p$add_argument('--celltype_label', type="character", help='Cell type label')
p$add_argument('--outdir',          type="character",                               help='Output file')

args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
args$metadata <- file.path(io$basedir,"results/mapping/sample_metadata_after_mapping.txt.gz")
args$outdir <- file.path(io$basedir,"results/celltype_proportions/fig")
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings = F)

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
  # "ExE_endoderm",
  "ExE_ectoderm"
  # "Parietal_endoderm"
)

opts$samples <- c(
  "WT_1", "WT_2", "WT_3", "WT_4", "WT_6",
  "Dnmt3a_KO_1", "Dnmt3a_KO_12", "Dnmt3a_KO_2", "Dnmt3a_KO_Dnmt3b_HET_1", "Dnmt3a_KO_Dnmt3b_HET_2",
  "Dnmt3b_KO_1", "Dnmt3b_KO_2", "Dnmt3b_KO_6", "Dnmt3b_KO_7", "Dnmt3b_KO_9", 
  "Dnmt1_KO_10", "Dnmt1_KO_15", "Dnmt1_KO_2", "Dnmt1_KO_3", "Dnmt1_KO_9"
  # "Dnmt3ab_KO_1", "Dnmt3ab_KO_2"
)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  # .[pass_rnaQC==TRUE & doublet_call==FALSE & !is.na(eval(as.name(args$celltype_label)))] %>%
  .[pass_rnaQC==TRUE & alias%in%opts$samples & celltype.mapped%in%opts$celltypes] %>%
  .[,alias:=factor(alias,levels=opts$samples)] %>%
  .[class=="Dnmt3a_KO_Dnmt3b_HET",class:="Dnmt3a_KO"] %>%
  setnames("celltype.mapped","celltype")

if (opts$remove_ExE_cells) {
  print("Removing ExE cells...")
  sample_metadata <- sample_metadata %>%
    .[!celltype%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

###############
## Per class ##
###############

to.plot <- sample_metadata %>%
  .[,N:=.N,by="class"] %>%
  .[,.(N=.N, celltype_proportion=.N/unique(N)),by=c("class","alias","celltype")] 

# Define colours and cell type order
opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors) %in% unique(to.plot$celltype)]
to.plot[,celltype:=factor(celltype, levels=rev(names(opts$celltype.colors)))]

# classes.to.plot <- unique(to.plot$class)

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
  
# if (length(unique(to.plot[class==i,alias]))>2) { width <- 12 } else { width <- 6 }

pdf(file.path(args$outdir,"celltype_proportions_horizontal_barplots.pdf"), width=10, height=14)
print(p)
dev.off()
