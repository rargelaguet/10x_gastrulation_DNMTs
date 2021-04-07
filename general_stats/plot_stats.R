#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
  source("/Users/ricard/10x_gastrulation_DNMTs/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/10x_gastrulation_DNMTs/settings.R")
  source("/homes/ricard/10x_gastrulation_DNMTs/utils.R")
}

io$outdir <- paste0(io$basedir,"/results/general_stats")

opts$samples <- c(
  # "E12.5_Dnmt3aHET_Dnmt3bWT_1",
  # "E12.5_Dnmt3aKO_Dnmt3bWT_1",
  # "E12.5_Dnmt3aWT_Dnmt3bHET_1",
  # "E12.5_Dnmt3aWT_Dnmt3bHET_2",
  # "E12.5_Dnmt3aWT_Dnmt3bKO_1",
  # "E12.5_Dnmt3aWT_Dnmt3bKO_2",
  "E8.5_Dnmt1KO_1",
  "E8.5_Dnmt1KO_2",
  "E8.5_Dnmt1KO_3",
  "E8.5_Dnmt3aHET_Dnmt3bKO_1",
  "E8.5_Dnmt3aHET_Dnmt3bWT_1",
  "E8.5_Dnmt3aKO_Dnmt3bHET_1",
  "E8.5_Dnmt3aKO_Dnmt3bHET_2",
  "E8.5_Dnmt3aKO_Dnmt3bKO_1",
  "E8.5_Dnmt3aKO_Dnmt3bKO_2",
  "E8.5_Dnmt3aKO_Dnmt3bWT_1",
  "E8.5_Dnmt3aKO_Dnmt3bWT_2",
  "E8.5_Dnmt3aWT_Dnmt3bKO_1",
  "E8.5_Dnmt3aWT_Dnmt3bKO_2",
  "E8.5_WT_1",
  "E8.5_WT_2",
  "E8.5_WT_3",
  "E8.5_WT_4",
  "E8.5_WT_5",
  "E8.5_WT_6",
  "E8.5_WT_7"
)

#########################
## Load sample metadata ##
#########################

sample_metadata <- fread(io$metadata) %>% 
  .[pass_QC==TRUE & sample%in%opts$samples]

###############################################
## Boxplots of general statistics per class ##
###############################################

to.plot <- sample_metadata %>% 
  melt(id.vars=c("cell","class"), measure.vars=c("nCount_RNA","nFeature_RNA"))

p <- ggboxplot(to.plot, x = "class", y = "value", outlier.shape=NA) +
  yscale("log10", .format = TRUE) +
  labs(x="", y="") +
  facet_wrap(~variable, scales="free", nrow=1) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    # axis.text.x = element_blank(),
    axis.text.x = element_text(colour="black",size=rel(0.9), angle=50, hjust=1),
    # axis.ticks.x = element_blank()
  )

pdf(paste0(io$outdir,"/general_stats_per_sample.pdf"), width=16, height=6, useDingbats = F)
print(p)
dev.off()

###########################################################
## Boxplots of general statistics per celltype and class ##
###########################################################

to.plot <- sample_metadata %>% 
  .[,.(N=.N, nFeature_RNA=mean(nFeature_RNA)),by=c("sample","class","celltype.mapped")] %>%
  .[N>15]

p <- ggboxplot(to.plot, x = "class", y = "nFeature_RNA", fill="class", outlier.shape=NA) +
  yscale("log10", .format = TRUE) +
  # scale_fill_manual(values=opts$celltype.colors) +
  labs(x="", y="") +
  facet_wrap(~celltype.mapped, scales="fixed") +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
    # axis.text.x = element_text(colour="black",size=rel(0.9), angle=50, hjust=1),
  )

pdf(paste0(io$outdir,"/general_stats_per_celltype_and_class.pdf"), width=10, height=10)
print(p)
dev.off()

########################################
## Barplots number of cells per sample ##
########################################

to.plot <- sample_metadata[,.N,by=c("sample","stage")]

p <- ggbarplot(to.plot, x = "sample", y = "N", fill="gray70") +
  labs(x="", y="Number of cells (after QC)") +
  facet_wrap(~stage, nrow=1, scales="free_x") +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    # axis.text.x = element_blank(),
    axis.text.x = element_text(colour="black",size=rel(0.8), angle=40, hjust=1),
    # axis.ticks.x = element_blank()
  )

pdf(paste0(io$outdir,"/N_per_sample.pdf"), width=8, height=6, useDingbats = F)
print(p)
dev.off()


##################################################
## Boxplots of general statistics per cell type ##
##################################################

to.plot <- sample_metadata %>% 
  melt(id.vars=c("cell","celltype.mapped"), measure.vars=c("nCount_RNA","nFeature_RNA"))

p <- ggboxplot(to.plot, x = "celltype.mapped", y = "value", fill="celltype.mapped", outlier.shape=NA) +
  yscale("log10", .format = TRUE) +
  labs(x="", y="") +
  scale_fill_manual(values=opts$celltype.colors) +
  facet_wrap(~variable, scales="free_y") +
  guides(fill = guide_legend(override.aes = list(size=0.25), ncol=1)) +
  scale_size(guide = 'none') +
  theme(
    legend.position = "right",
    legend.text = element_text(size=rel(0.75)),
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

pdf(paste0(io$outdir,"/general_stats_per_celltype.pdf"), width=16, height=10, useDingbats = F)
print(p)
dev.off()

