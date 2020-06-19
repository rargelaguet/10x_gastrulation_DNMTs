library(ggpubr)

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/10x_gastrulation_DNMTs/settings.R")
}

io$outdir <- paste0(io$basedir,"/results/rna_stats")

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
  
  "E125_DNMT3A_HET_A_L001",
  "E125_DNMT3A_HET_A_L003",
  "E125_DNMT3A_KO_B_L002",
  "E125_DNMT3A_KO_E_L004",
  "A_E12_5_D3a_Het_L001",
  "B_E12_5_D3a_KO_L002"
)

##########
## Plot ##
##########

# sample_metadata <- fread(io$metadata)

to.plot <- sample_metadata %>% 
  .[batch%in%opts$batches] %>%
  melt(id.vars=c("cell","batch","class","stage"), measure.vars=c("nCount_RNA","nFeature_RNA"))

p <- ggboxplot(to.plot, x = "batch", y = "value", fill="stage", outlier.shape=NA) +
  yscale("log10", .format = TRUE) +
  labs(x="", y="") +
  facet_wrap(~variable, scales="free_y") +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# pdf(paste0(io$outdir,"/rna_stats.pdf"), width=6, height=4, useDingbats = F)
print(p)
# dev.off()
