#####################
## Define settings ##
#####################

source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
# source("/Users/ricard/10x_gastrulation_DNMTs/mapping/plot/plot_utils.R")

################
## Define I/O ##
################

io$outdir <- paste0(io$basedir,"/results/celltype_proportions/barplots")
dir.create(paste0(io$outdir,"/per_embryo"), showWarnings = F)
dir.create(paste0(io$outdir,"/per_class"), showWarnings = F)

####################
## Define options ##
####################

opts$batches <- c(
  "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001",
  "SIGAB6_E85_3_Dnmt3aWT_Dnmt3b_WT_L002",
  "SIGAC6_E85_5_Dnmt3aKO_Dnmt3b_Het_L003",
  "SIGAD6_E85_8_Dnmt3aHet_Dnmt3b_KO_L004",
  "15_E8_5_D3A_WT_D3B_WT_L007",
  "17_E8_5_D3A_KO_D3B_WT_L008",
  "2_E8_5_D3A_WT_D3B_KO_L003",
  "3_E8_5_D3A_HET_D3B_WT_L004",
  "7_E8_5_D3A_WT_D3B_KO_L005",
  "8_E8_5_D3A_KO_D3B_KO_L006",
  "E8_5_Dnmt1_KO_male_SIGAC8_L001",
  "E8_5_Dnmt1_KO_male_SIGAD8_L002",
  "E8_5_Dnmt1_KO_male_SIGAE8_L003",
  "E8_5_Dnmt1_WT_female_SIGAB8_L004",
  "E8_5_Dnmt1_WT_female_SIGAF8_L005",
  "E8_5_Dnmt3ab_WT_female_SIGAA8_L006",
  "SIGAH10_Dnmt3ab_WT_L002",
  "SIGAH11_Dnmt3ab_WT_L003",
  "SIGAH9_Dnmt3a_KO_Dnmt3b_Het_L001",
  "SIGAG5_9_dnmt3ab_DKO_L005"
)

###############
## Load data ##
###############

# mapping <- fread(io$mapping)
# sample_metadata <- sample_metadata %>% merge(mapping)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  .[pass_QC==TRUE & !is.na(celltype.mapped)] %>%
  .[batch%in%opts$batches]
  
# Rename samples
foo <- sample_metadata[,c("batch","class")] %>% unique %>% .[,sample:=paste(class,1:.N,sep="_"),by="class"]
sample_metadata <- sample_metadata %>% merge(foo,by=c("batch","class"))

# Define cell type order
opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors) %in% sample_metadata$celltype.mapped]
stopifnot(sort(unique(as.character(sample_metadata$celltype.mapped))) == sort(names(opts$celltype.colors)))
sample_metadata <- sample_metadata %>% 
  .[,celltype.mapped:=factor(celltype.mapped,levels=sort(names(opts$celltype.colors), decreasing = F))]

#####################################
## Calculate cell type proportions ##
#####################################

to.plot.embryo <- sample_metadata %>% 
  .[,.N, by=c("celltype.mapped","sample","class")]# %>%
  # .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped," ", "_")] %>

to.plot.class <- sample_metadata %>%
  .[,.N, by=c("celltype.mapped","class")]


#####################
## Plot per class ##
#####################

for (i in unique(to.plot.class$class)) {
  
  p <- ggplot(to.plot.class[class==i], aes(x=celltype.mapped, y=N)) +
    geom_bar(aes(fill=celltype.mapped), stat="identity", color="black") +
    scale_fill_manual(values=opts$celltype.colors, drop=FALSE) +
    scale_x_discrete(drop=FALSE) +
    coord_flip() +
    labs(y="Number of cells") +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=rel(1.2)),
      axis.title.x = element_text(color="black", size=rel(1.1)),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1.1), color="black"),
      axis.text.x = element_text(size=rel(1.1), color="black")
    )
  
  pdf(sprintf("%s/per_class/barplots_%s.pdf",io$outdir,i), width=7, height=7)
  print(p)
  dev.off()
  
  p <- ggplot(to.plot.embryo[class==i], aes(x=celltype.mapped, y=N)) +
    geom_bar(aes(fill=celltype.mapped), stat="identity", color="black") +
    scale_fill_manual(values=opts$celltype.colors, drop=FALSE) +
    scale_x_discrete(drop=FALSE) +
    facet_wrap(~sample, nrow=1) +
    coord_flip() +
    labs(y="Number of cells") +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=rel(0.9)),
      axis.title.x = element_text(color="black", size=rel(1.1)),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1.1), color="black"),
      axis.text.x = element_text(size=rel(1.1), color="black")
    )
  
  pdf(sprintf("%s/per_class/barplots_facet_embryo_%s.pdf",io$outdir,i))
  print(p)
  dev.off()
}

#####################
## Plot per embryo ##
#####################

for (i in unique(to.plot.embryo$sample)) {
  
  p <- ggplot(to.plot.embryo[sample==i], aes(x=celltype.mapped, y=N)) +
    geom_bar(aes(fill=celltype.mapped), stat="identity", color="black") +
    scale_fill_manual(values=opts$celltype.colors, drop=FALSE) +
    scale_x_discrete(drop=FALSE) +
    coord_flip() +
    labs(y="Number of cells") +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=rel(1.2)),
      axis.title.x = element_text(color="black", size=rel(1.1)),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1.1), color="black"),
      axis.text.x = element_text(size=rel(1.1), color="black")
    )
  
  pdf(sprintf("%s/per_embryo/barplots_%s.pdf",io$outdir,i), width=7, height=7)
  print(p)
  dev.off()
}
