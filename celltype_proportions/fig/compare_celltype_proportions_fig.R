here::i_am("celltype_proportions/compare_celltype_proportions.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

io$metadata <- file.path(io$basedir,"results/mapping/sample_metadata_after_mapping.txt.gz")
io$outdir <- file.path(io$basedir,"results/celltype_proportions/comparisons/fig"); dir.create(io$outdir, showWarnings = F)

####################
## Define options ##
####################

opts$ko.classes <- c(
  "Dnmt3a_KO", 
  "Dnmt3b_KO",
  "Dnmt1_KO"
  # "Dnmt3ab_KO"
)

opts$wt.classes <- c("WT")
  
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

opts$remove.ExE.celltypes <- FALSE
opts$remove.blood <- FALSE
opts$remove.small.lineages <- FALSE
opts$remove.small.embryos <- TRUE

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE & celltype.mapped%in%opts$celltypes & class%in%c(opts$wt.classes,opts$ko.classes)] %>%
  .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,opts$rename_celltypes)] %>%
  .[,celltype.mapped:=factor(celltype.mapped,levels=unique(celltype.mapped))] %>%
  .[,dataset:=ifelse(grepl("Grosswendt",sample),"CRISPR","KO")] %>%
  .[,c("cell","sample","alias","class","celltype.mapped","dataset")]

# Filter cells
if (opts$remove.blood) {
  sample_metadata <- sample_metadata %>% .[!celltype.mapped=="Erythroid"]
}

if (opts$remove.ExE.celltypes) {
  sample_metadata <- sample_metadata %>%
    # .[!celltype.mapped%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")] %>%
    .[!celltype.mapped%in%c("ExE_ectoderm","Parietal_endoderm")]
}

# if (opts$remove.small.lineages) {
#   opts$min.cells <- 250
#   sample_metadata <- sample_metadata %>%
#     .[,N:=.N,by=c("celltype.mapped")] %>% .[N>opts$min.cells] %>% .[,N:=NULL]
# }

if (opts$remove.small.embryos) {
  opts$min.cells <- 1000
  sample_metadata <- sample_metadata %>%
    .[,N:=.N,by="alias"] %>% .[N>opts$min.cells] %>% .[,N:=NULL]
}

# Print statistics
print(table(sample_metadata$sample))
print(table(sample_metadata$celltype.mapped))

####################################
## Calculate celltype proportions ##
####################################

# Calculate celltype proportions for WT
wt_proportions.dt <- sample_metadata %>%
  .[class%in%opts$wt.classes] %>%
  .[,ncells:=.N] %>%
  .[,.(proportion=.N/unique(ncells), N=round(.N/length(unique(sample)))), by="celltype.mapped"]


# Calculate celltype proportions for KO samples
ko_proportions_per_sample.dt <- sample_metadata %>%
  .[class%in%opts$ko.classes] %>%
  setkey(celltype.mapped,sample) %>%
  .[CJ(celltype.mapped,sample, unique = TRUE), .N, by = .EACHI] %>%
  merge(unique(sample_metadata[,c("sample","class")]), by="sample") %>%
  .[,ncells:=sum(N), by="sample"] %>% .[,proportion:=(N+1)/ncells]

# Calculate celltype proportions for KO classes
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

# Define order based on Dnmt1 KO
celltype.order <- proportions_per_class.dt[class=="Dnmt1_KO"] %>% setorder(-diff_proportion) %>% .$celltype.mapped %>% as.character

# For viz purposes
proportions_per_sample.dt[diff_proportion>=6.5,diff_proportion:=6.5]
# ylimits <- max(abs(proportions_per_sample.dt$diff_proportion)) + 0.25
ylimits <- 6.75

# i <- "Dnmt1_KO"
for (i in opts$ko.classes) {
  
  celltypes.to.plot <- proportions_per_sample.dt %>%
    .[class==i,.(N=sum(N.ko)+sum(N.wt)),by=c("class","celltype.mapped")] %>% 
    .[N>=50,celltype.mapped] %>% as.character
  
  to.plot <- proportions_per_sample.dt %>%
    .[class==i & celltype.mapped%in%celltypes.to.plot] %>%
    merge(unique(sample_metadata[,c("sample","dataset")]),by="sample")
  
  # celltype.order <- to.plot %>% .[,mean(diff_proportion),by="celltype.mapped"] %>% setorder(-V1) %>% .$celltype.mapped
  to.plot <- to.plot %>% .[,celltype.mapped:=factor(celltype.mapped,levels=celltype.order)]
  
  # text.dt <- proportions_per_class.dt %>% 
  #   .[class==i & celltype.mapped%in%celltype.order] %>% 
  #   .[,celltype.mapped:=factor(celltype.mapped,levels=celltype.order)]
  
  p <- ggplot(to.plot, aes(x=celltype.mapped, y=diff_proportion)) +
    geom_point(aes(fill = celltype.mapped), shape=21, size=1.25, stroke=0.1) +
    geom_boxplot(aes(fill = celltype.mapped), alpha=0.75) +
    facet_wrap(~dataset) +
    # geom_text(y=-ylimits, aes(label=N.wt), size=2.5, data=text.dt) +
    # geom_text(y=ylimits, aes(label=N.ko), size=2.5, data=text.dt) +
    coord_flip(ylim=c(-ylimits,ylimits)) +
    geom_hline(yintercept=0, linetype="dashed", size=0.5) +
    scale_fill_manual(values=opts$celltype.colors, drop=F) +
    theme_classic() +
    labs(y="Difference in proportions (log2)", x="", title=i) +
    theme(
      legend.position = "none",
      # axis.title = element_blank(),
      plot.title = element_text(size=rel(1.25), hjust=0.5, color="black"),
      axis.text.y = element_text(color="black"),
      axis.text.x = element_text(color="black")
    )
  
  # pdf(file.path(io$outdir,sprintf("%s_boxplots.pdf",i)), width=6, height=7)
  pdf(file.path(io$outdir,sprintf("%s_boxplots_per_dataset.pdf",i)), width=9, height=5)
  print(p)
  dev.off()
}

########################################
## Plot correlation between data sets ##  
########################################

wt_proportions_dataset.dt <- sample_metadata %>%
  .[class%in%opts$wt.classes] %>%
  .[,ncells:=.N,by="dataset"] %>%
  .[,.(proportion=.N/unique(ncells), N=round(.N/length(unique(sample)))), by=c("celltype.mapped","dataset")]

ko_proportions_per_class_dataset.dt <- sample_metadata %>%
  .[!class%in%opts$wt.classes] %>%
  setkey(celltype.mapped,class,dataset) %>%
  .[CJ(celltype.mapped,class,dataset,unique = TRUE), .N, by = .EACHI] %>%
  .[,ncells:=sum(N), by=c("class","dataset")] %>% .[,proportion:=(N+1)/ncells]

to.plot <- merge(
  ko_proportions_per_class_dataset.dt, 
  wt_proportions_dataset.dt, 
  by = c("celltype.mapped","dataset"), allow.cartesian=T, suffixes = c(".ko",".wt")) %>% 
  .[,N:=N.ko+N.wt] %>% .[N>=50] %>%
  .[,c("diff_proportion"):=list(log2(proportion.ko/proportion.wt))] %>%
  dcast(celltype.mapped+class~dataset,value.var="diff_proportion")


p <- ggscatter(to.plot, x="CRISPR", y="KO", fill="celltype.mapped", shape=21, size=4,
               add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
  scale_fill_manual(values=opts$celltype.colors) +
  stat_cor(method = "pearson") +
  facet_wrap(~class, scales="fixed") +
  theme_classic() +
  labs(y="log2 difference in proportions (CRISPR)", x="log2 difference in proportions (KO)") +
  theme(
    legend.position = "none",
    axis.text = element_text(color="black", size=rel(0.90))
  )

pdf(file.path(io$outdir,"scatterplot_diff_proportions.pdf"), width=8.5, height=5)
print(p)
dev.off()
