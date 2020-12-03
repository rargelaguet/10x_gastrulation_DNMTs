linMap <- function(x, from, to) return( (x - min(x)) / max(x - min(x)) * (to - from) + from )

source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")

################
## Define I/O ##
################

io$outdir <- paste0(io$basedir,"/results/mapping/pdf")

####################
## Define options ##
####################

opts$classes <- c(
  # "E12.5_Dnmt3aWT_Dnmt3bHET",
  # "E12.5_Dnmt3aWT_Dnmt3bKO",
  # "E12.5_Dnmt3aHET_Dnmt3bWT",
  # "E12.5_Dnmt3aKO_Dnmt3bWT",
  "E8.5_Dnmt3aKO_Dnmt3bWT",
  "E8.5_WT",
  "E8.5_Dnmt3aHET_Dnmt3bKO",
  "E8.5_Dnmt3aHET_Dnmt3bWT",
  "E8.5_Dnmt3aKO_Dnmt3bHET",
  "E8.5_Dnmt3aKO_Dnmt3bKO",
  "E8.5_Dnmt3aWT_Dnmt3bKO",
  "E8.5_Dnmt1KO"
  )

opts$to.merge <- c(
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

############################
## Update sample metadata ##
############################

sample_metadata <- fread(io$metadata) %>%
  .[pass_QC==TRUE] %>%
  .[class%in%opts$classes & !is.na(celltype.mapped)] %>%
  .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,opts$to.merge)]# %>%
  # .[!celltype.mapped%in%c("Erythroid")]

table(sample_metadata$batch)


################
## Parse data ##
################

matrix <- sample_metadata %>%
  .[,ncells:=.N,by="batch"] %>%
  .[,.(proportion=.N/unique(ncells), N=.N),by=c("celltype.mapped","batch")] %>%
  dcast(batch~celltype.mapped, fill=0, value.var="proportion") %>%
  matrix.please

#########
## PCA ##
#########

pca <- prcomp(matrix, rank.=5)

variance.explained.by.pc <- 100*(pca$sdev / sum(pca$sdev))

plot(variance.explained.by.pc)

##########
## Plot ##
##########

to.plot <- pca$rotation %>% as.data.table %>% 
  .[,celltype:=rownames(pca$rotation)]

p <- ggplot(to.plot, aes(x=PC1, y=PC3, fill=celltype)) +
  geom_point(shape=21, stroke=0.5, color="black", size=5) +
  scale_fill_manual(values=opts$celltype.colors) +
  # labs(x=sprintf("PC1 (%.2f%%)",variance.explained.by.pc[1]), y=sprintf("PC2 (%.2f%%)",variance.explained.by.pc[2])) +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    # axis.text = element_blank(),
    axis.text = element_text(color="black", size=rel(0.75)),
    axis.ticks = element_blank()
  )

# pdf(paste0(io$outdir,"/pca_mapping_stages.pdf"), width=7, height=5, useDingbats = F)
print(p)
# dev.off()

##########
## Plot ##
##########

to.plot <- pca$x %>% as.data.table %>% 
  .[,batch:=rownames(pca$x)] %>%
  merge(unique(sample_metadata[,c("batch","class")]), by="batch")

p <- ggplot(to.plot, aes(x=PC1, y=PC4, fill=class)) +
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
