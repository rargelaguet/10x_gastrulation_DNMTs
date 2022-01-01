 here::i_am("celltype_proportions/pca_celltype_proportions.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

################
## Define I/O ##
################

io$outdir <- paste0(io$basedir,"/results_new/celltype_proportions")

####################
## Define options ##
####################

opts$classes <- c(
  "E8.5_Dnmt3aKO_Dnmt3bWT",
  "E8.5_WT",
  "E8.5_Dnmt3aHET_Dnmt3bKO",
  # "E8.5_Dnmt3aHET_Dnmt3bWT",  # THIS ONE WAS COMMENTED OUT
  "E8.5_Dnmt3aKO_Dnmt3bHET",
  "E8.5_Dnmt3aKO_Dnmt3bKO",
  "E8.5_Dnmt3aWT_Dnmt3bKO",
  "E8.5_Dnmt1KO"
  )

opts$rename_celltypes <- c(
  "Erythroid3" = "Erythroid",
  "Erythroid2" = "Erythroid",
  "Erythroid1" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors",
  "Intermediate_mesoderm" = "Mixed_mesoderm",
  "Paraxial_mesoderm" = "Mixed_mesoderm",
  "Nascent_mesoderm" = "Mixed_mesoderm",
  "Pharyngeal_mesoderm" = "Mixed_mesoderm",
  "Visceral_endoderm" = "ExE_endoderm"
)

opts$remove.ExE.celltypes <- TRUE
opts$remove.blood <- TRUE
opts$remove.small.lineages <- FALSE

############################
## Update sample metadata ##
############################

sample_metadata <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE & class%in%opts$classes & !is.na(celltype.mapped)] %>%
  .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,opts$rename_celltypes)]

# Filter cells
if (opts$remove.blood) {
  sample_metadata <- sample_metadata %>% .[!celltype.mapped=="Erythroid"]
}
if (opts$remove.ExE.celltypes) {
  sample_metadata <- sample_metadata %>%
    .[!celltype.mapped%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
    # .[!celltype.mapped%in%c("ExE_ectoderm","Parietal_endoderm")]
}
if (opts$remove.small.lineages) {
  opts$min.cells <- 100
  sample_metadata <- sample_metadata %>%
    .[,N:=.N,by=c("celltype.mapped")] %>% .[N>opts$min.cells] %>% .[,N:=NULL]
}

# print stats
table(sample_metadata$alias)

################
## Parse data ##
################

mtx <- sample_metadata %>%
  .[,ncells:=.N,by="alias"] %>%
  .[,.(proportion=.N/unique(ncells), N=.N),by=c("celltype.mapped","alias")] %>%
  dcast(alias~celltype.mapped, fill=0, value.var="proportion") %>%
  matrix.please

#########
## PCA ##
#########

pca <- prcomp(mtx, rank.=5)

variance.explained.by.pc <- 100*(pca$sdev / sum(pca$sdev))

plot(variance.explained.by.pc)

##################
## Plot samples ##
##################

to.plot <- pca$x %>% as.data.table %>% 
  .[,alias:=rownames(pca$x)] %>%
  merge(unique(sample_metadata[,c("alias","class")]), by="alias")

p <- ggplot(to.plot, aes(x=PC1, y=PC2, fill=class)) +
  geom_point(shape=21, stroke=0.5, color="black", size=5) +
  scale_fill_brewer(palette="Dark2") +
  # scale_shape_manual(values=c(21,24)) +
  # guides(fill=guide_legend(override.aes=list(shape=21))) +
  # guides(shape=guide_legend(override.aes=list(fill="black"))) +
  # labs(x=sprintf("PC1 (%.2f%%)",variance.explained.by.pc[1]), y=sprintf("PC2 (%.2f%%)",variance.explained.by.pc[2])) +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    axis.text = element_blank(),
    # axis.text = element_text(color="black", size=rel(0.8))
    axis.ticks = element_blank(),
  )

# pdf(paste0(io$outdir,"/pca_mapping_stages.pdf"), width=7, height=5, useDingbats = F)
print(p)
# dev.off()

##########################
## Plot feature weights ##
##########################

to.plot <- pca$rotation %>% as.data.table %>% 
  .[,celltype:=rownames(pca$rotation)]

p <- ggplot(to.plot, aes(x=PC1, y=PC2, fill=celltype)) +
  geom_point(shape=21, stroke=0.5, color="black", size=5) +
  scale_fill_manual(values=opts$celltype.colors) +
  # labs(x=sprintf("PC1 (%.2f%%)",variance.explained.by.pc[1]), y=sprintf("PC2 (%.2f%%)",variance.explained.by.pc[2])) +
  theme_classic() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    # axis.text = element_blank(),
    axis.text = element_text(color="black", size=rel(0.75)),
    axis.ticks = element_blank()
  )

# pdf(paste0(io$outdir,"/pca_mapping_stages.pdf"), width=7, height=5, useDingbats = F)
print(p)
# dev.off()

