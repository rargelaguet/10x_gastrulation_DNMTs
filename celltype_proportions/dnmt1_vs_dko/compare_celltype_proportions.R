here::i_am("celltype_proportions/compare_celltype_proportions.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--celltype_label', type="character", help='Cell type label')
p$add_argument('--outdir',          type="character",                               help='Output file')

args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
args$metadata <- file.path(io$basedir,"results_new/mapping/sample_metadata_after_mapping.txt.gz")
args$celltype_label <- "celltype.mapped"
args$outdir <- file.path(io$basedir,"results_new/celltype_proportions/comparisons/dnmt1_vs_dnmt3ab_ko")
## END TEST ##

# I/O
dir.create(file.path(args$outdir,"boxplots"), showWarnings = F)
dir.create(file.path(args$outdir,"boxplots/per_class"), showWarnings = F)
dir.create(file.path(args$outdir,"boxplots/per_sample"), showWarnings = F)
dir.create(file.path(args$outdir,"polar_plots"), showWarnings = F)
dir.create(file.path(args$outdir,"polar_plots/per_class"), showWarnings = F)
dir.create(file.path(args$outdir,"polar_plots/per_sample"), showWarnings = F)

####################
## Define options ##
####################

opts$classes <- c(
  "E8.5_Dnmt1KO",
  "E8.5_Dnmt3aKO_Dnmt3bKO"
)

opts$wt.classes <- c("E8.5_Dnmt1KO")
  
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
  "ExE_ectoderm"
  # "Parietal_endoderm"
)

opts$rename_celltypes <- c(
  "Erythroid3" = "Erythroid",
  "Erythroid2" = "Erythroid",
  "Erythroid1" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors"
  # "Intermediate_mesoderm" = "Mixed_mesoderm",
  # "Paraxial_mesoderm" = "Mixed_mesoderm",
  # "Nascent_mesoderm" = "Mixed_mesoderm",
  # "Pharyngeal_mesoderm" = "Mixed_mesoderm"
  # "Visceral_endoderm" = "ExE_endoderm"
)

opts$remove.ExE.celltypes <- TRUE
opts$remove.blood <- FALSE
opts$remove.small.lineages <- FALSE

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_rnaQC==TRUE & celltype.mapped%in%opts$celltypes & class%in%opts$classes] %>%
  .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,opts$rename_celltypes)] %>%
  .[,celltype.mapped:=factor(celltype.mapped,levels=unique(celltype.mapped))] %>%
  .[,c("cell","sample","alias","class","celltype.mapped")]

# Filter cells
if (opts$remove.blood) {
  sample_metadata <- sample_metadata %>% .[!celltype.mapped=="Erythroid"]
}

if (opts$remove.ExE.celltypes) {
  sample_metadata <- sample_metadata %>%
    # .[!celltype.mapped%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")] %>%
    .[!celltype.mapped%in%c("ExE_ectoderm","Parietal_endoderm")]
}

if (opts$remove.small.lineages) {
  opts$min.cells <- 250
  sample_metadata <- sample_metadata %>%
    .[,N:=.N,by=c("celltype.mapped")] %>% .[N>opts$min.cells] %>% .[,N:=NULL]
}

# Print statistics
print(table(sample_metadata$class))
print(table(sample_metadata$celltype.mapped))

####################################
## Calculate celltype proportions ##
####################################

# Calculate celltype proportions for WT
# wt_proportions.dt <- sample_metadata %>%
#   .[class%in%opts$wt.classes] %>%
#   setkey(celltype.mapped,class) %>%
#   .[CJ(celltype.mapped,class, unique = TRUE), .N, by = .EACHI] %>%
#   .[,ncells:=round(mean(N)), by="class"] %>% .[,proportion:=(N+1)/ncells]
wt_proportions.dt <- sample_metadata %>%
  .[class%in%opts$wt.classes] %>%
  .[,ncells:=.N] %>%
  .[,.(proportion=.N/unique(ncells), N=round(.N/length(unique(sample)))), by="celltype.mapped"]


# Calculate celltype proportions for KO samples
# ko_proportions_per_sample.dt <- sample_metadata %>%
#   .[!class%in%opts$wt.classes] %>%
#   .[,ncells:=.N, by="sample"] %>%
#   .[,.(proportion=.N/unique(ncells), N=.N),by=c("celltype.mapped","sample","class")]
ko_proportions_per_sample.dt <- sample_metadata %>%
  .[!class%in%opts$wt.classes] %>%
  setkey(celltype.mapped,sample) %>%
  .[CJ(celltype.mapped,sample, unique = TRUE), .N, by = .EACHI] %>%
  merge(unique(sample_metadata[,c("sample","class")]), by="sample") %>%
  .[,ncells:=sum(N), by="sample"] %>% .[,proportion:=(N+1)/ncells]

# Calculate celltype proportions for KO classes
# ko_proportions_per_class.dt <- sample_metadata %>%
#   .[!class%in%opts$wt.classes] %>%
#   .[,ncells:=.N, by="class"] %>%
#   .[,.(proportion=.N/unique(ncells), N=.N),by=c("celltype.mapped","class")]
ko_proportions_per_class.dt <- sample_metadata %>%
  .[!class%in%opts$wt.classes] %>%
  setkey(celltype.mapped,class) %>%
  .[CJ(celltype.mapped,class, unique = TRUE), .N, by = .EACHI] %>%
  .[,ncells:=sum(N), by="class"] %>% .[,proportion:=(N+1)/ncells]

# Merge
proportions_per_sample.dt <- merge(
  ko_proportions_per_sample.dt, 
  wt_proportions.dt, 
  by = c("celltype.mapped"), allow.cartesian=T, suffixes = c(".ko",".wt")
) %>% .[,c("diff_proportion"):=list(log2(proportion.ko/proportion.wt))] 

proportions_per_class.dt <- merge(
  ko_proportions_per_class.dt, 
  wt_proportions.dt, 
  by = c("celltype.mapped"), allow.cartesian=T, suffixes = c(".ko",".wt")
) %>% .[,c("diff_proportion"):=list(log2(proportion.ko/proportion.wt))] 


########################
## Boxplots per class ##  
#######################

ylimits <- max(abs(proportions_per_sample.dt$diff_proportion)) + 0.25

for (i in unique(proportions_per_sample.dt$class)) {
  
  celltypes.to.plot <- proportions_per_sample.dt %>%
    .[class==i,.(N=sum(N.ko)+sum(N.wt)),by=c("class","celltype.mapped")] %>% 
    .[N>=50,celltype.mapped] %>% as.character
  
  to.plot <- proportions_per_sample.dt %>%
    .[class==i & celltype.mapped%in%celltypes.to.plot] 
  
  celltype.order <- to.plot %>%
    .[,mean(diff_proportion),by="celltype.mapped"] %>% setorder(-V1) %>% .$celltype.mapped
  to.plot <- to.plot %>% .[,celltype.mapped:=factor(celltype.mapped,levels=celltype.order)]
  
  # text.dt <- proportions_per_class.dt %>% 
  #   .[class==i & celltype.mapped%in%celltype.order] %>% 
  #   .[,celltype.mapped:=factor(celltype.mapped,levels=celltype.order)]
  
  p <- ggplot(to.plot, aes(x=celltype.mapped, y=diff_proportion)) +
    geom_point(aes(fill = celltype.mapped), shape=21, size=1) +
    geom_boxplot(aes(fill = celltype.mapped), alpha=0.5) +
    # geom_text(y=-ylimits, aes(label=N.wt), size=2.5, data=text.dt) +
    # geom_text(y=ylimits, aes(label=N.ko), size=2.5, data=text.dt) +
    coord_flip(ylim=c(-ylimits,ylimits)) +
    geom_hline(yintercept=0, linetype="dashed", size=0.5) +
    scale_fill_manual(values=opts$celltype.colors, drop=F) +
    theme_classic() +
    labs(y="Difference in proportions (log2)", x="") +
    theme(
      legend.position = "none",
      # axis.title = element_blank(),
      axis.text.y = element_text(color="black"),
      axis.text.x = element_text(color="black")
    )
  
  # pdf(sprintf("%s/boxplots/per_class/%s_boxplots.pdf",args$outdir,i))
  print(p)
  # dev.off()
}


###########################
## Polar plots per class ##
###########################

# ylimits <- max(abs(proportions_per_class.dt[!is.infinite(diff_proportion),diff_proportion]))
ylimits <- 6

to.plot.wt_line <- data.table(
  celltype.mapped = unique(proportions_per_sample.dt$celltype.mapped),
  diff_proportion = log2(1)
)

for (i in unique(proportions_per_class.dt$class)) {
  
  celltypes.to.plot <- proportions_per_sample.dt %>%
    .[class==i,.(N=sum(N.ko)+sum(N.wt)),by=c("class","celltype.mapped")] %>% 
    .[N>=25,celltype.mapped] %>% as.character
  
  to.plot <- proportions_per_sample.dt %>%
    .[class==i & celltype.mapped%in%celltypes.to.plot]  %>%
    .[diff_proportion>=ylimits,diff_proportion:=ylimits] %>%
    .[diff_proportion<=(-ylimits),diff_proportion:=(-ylimits)]
  
  p <- ggplot(to.plot, aes(x=factor(celltype.mapped), y=-diff_proportion, group=1)) +
    geom_jitter(aes(fill = celltype.mapped), size=3.5, stat = 'identity', shape=21, width=0.1) +
    geom_polygon(color="black", fill=NA, alpha=0.5, linetype="dashed", data=to.plot.wt_line[celltype.mapped%in%to.plot$celltype.mapped]) +
    scale_fill_manual(values=opts$celltype.colors, drop=F) +
    coord_polar() + ylim(ylimits+0.1,-ylimits-0.1) +
    guides(colour = guide_legend(override.aes = list(size=2), ncol=1)) +
    scale_size(guide = 'none') +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.text =  element_text(size=rel(0.5)),
      legend.text = element_text(size=rel(0.75)),
      legend.title = element_blank(),
      axis.title=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      axis.line=element_blank(),
      axis.text.x = element_blank()
    )
  
  # pdf(sprintf("%s/polar_plots/per_class/%s_polar_plot.pdf",args$outdir,i), width = 5, height=5)
  print(p)
  # dev.off()
}
