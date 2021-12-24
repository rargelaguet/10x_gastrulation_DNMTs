#####################
## Define settings ##
#####################

source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
# source("/Users/ricard/10x_gastrulation_DNMTs/mapping/plot/plot_utils.R")

################
## Define I/O ##
################

io$metadata <- "/Users/ricard/data/10x_gastrulation_DNMTs/processed/sixth_batch/metadata.txt.gz"
io$outdir <- paste0(io$basedir,"/results/mapping/pdf")

####################
## Define options ##
####################

opts$batches <- c(
  # E12.5  
  # "E125_DNMT3A_HET_A_L001",
  # "E125_DNMT3A_HET_A_L003",
  # "E125_DNMT3A_KO_B_L002",
  # "E125_DNMT3A_KO_E_L004",
  # "A_E12_5_D3a_Het_L001",
  # "B_E12_5_D3a_KO_L002",
  
  
  # E8.5  
  # "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001",
  # "SIGAB6_E85_3_Dnmt3aWT_Dnmt3b_WT_L002",
  # "SIGAC6_E85_5_Dnmt3aKO_Dnmt3b_Het_L003",
  # "SIGAD6_E85_8_Dnmt3aHet_Dnmt3b_KO_L004",
  # "15_E8_5_D3A_WT_D3B_WT_L007",
  # "17_E8_5_D3A_KO_D3B_WT_L008",
  # "2_E8_5_D3A_WT_D3B_KO_L003",
  # "3_E8_5_D3A_HET_D3B_WT_L004",
  # "7_E8_5_D3A_WT_D3B_KO_L005",
  # "8_E8_5_D3A_KO_D3B_KO_L006",
  # "E8_5_Dnmt1_KO_male_SIGAC8_L001",
  # "E8_5_Dnmt1_KO_male_SIGAD8_L002",
  # "E8_5_Dnmt1_KO_male_SIGAE8_L003",
  # "E8_5_Dnmt1_WT_female_SIGAB8_L004",
  # "E8_5_Dnmt1_WT_female_SIGAF8_L005",
  # "E8_5_Dnmt3ab_WT_female_SIGAA8_L006",
  # "SIGAH10_Dnmt3ab_WT_L002",
  # "SIGAH11_Dnmt3ab_WT_L003",
  # "SIGAH9_Dnmt3a_KO_Dnmt3b_Het_L001"
  "SIGAG5_9_dnmt3ab_DKO_L005"
)

###############
## Load data ##
###############

# mapping <- fread(io$mapping)
# sample_metadata <- sample_metadata %>% merge(mapping)

################
## Parse data ##
################

to.plot <- fread(io$metadata) %>% 
  .[pass_QC==TRUE] %>%
  # .[class%in%opts$classes] %>%
  .[batch%in%opts$batches] %>%
  .[!is.na(celltype.mapped),.N, by=c("stage","celltype.mapped","batch","class")] %>%
  # .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped," ", "_")] %>%
  .[,celltype.mapped:=factor(celltype.mapped,levels=opts$celltypes)]


##########
## Plot ##
##########

for (i in unique(to.plot$batch)) {
  p <- ggplot(to.plot[batch==i], aes(x=celltype.mapped, y=N)) +
    geom_bar(aes(fill=celltype.mapped), stat="identity", color="black") +
    scale_fill_manual(values=opts$celltype.colors, drop=FALSE) +
    scale_x_discrete(drop=FALSE) +
    # facet_wrap(~batch, nrow=1, scales="fixed") +
    coord_flip() +
    labs(y="Number of cells") +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=rel(1.2)),
      axis.title.x = element_text(color="black", size=rel(1.1)),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1.3), color="black"),
      axis.text.x = element_text(size=rel(1.1), color="black")
    )
  
  pdf(sprintf("%s/barplots_%s.pdf",io$outdir,i), width=9, height=7)
  print(p)
  dev.off()
}
