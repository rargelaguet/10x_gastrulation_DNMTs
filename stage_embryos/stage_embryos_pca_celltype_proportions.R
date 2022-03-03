here::i_am("stage_embryos/pca_celltype_proportions.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$outdir <- file.path(io$basedir,"results/stage_embryos"); dir.create(io$outdir, showWarnings = F)

# Options
opts$classes <- c(
  # "E12.5_Dnmt3aWT_Dnmt3bHET",
  # "E12.5_Dnmt3b_KO",
  # "E12.5_Dnmt3a_HET_Dnmt3b_WT",
  # "E12.5_Dnmt3a_KO",
  # "Dnmt3a_HET_Dnmt3b_KO",
  # "Dnmt3a_HET_Dnmt3b_WT",
  # "Dnmt3a_KO_Dnmt3b_HET",
  "WT",
  "Dnmt3a_KO",
  "Dnmt3b_KO",
  "Dnmt1_KO"
  # "Dnmt3ab_KO"
)

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

opts$remove.ExE.celltypes <- FALSE

##########################
## Load sample metadata ##
##########################

cell_metadata.dt <- fread(io$metadata) %>% 
  .[,dataset:=ifelse(grepl("Grosswendt",sample),"CRISPR","KO")] %>% .[,dataset:=factor(dataset,levels=c("KO","CRISPR"))] %>%
  .[pass_rnaQC==TRUE & class%in%opts$classes]# %>%
  # .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped," ","_")] %>%
  # .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,"/","_")]

if (opts$remove.ExE.celltypes) {
  cell_metadata.dt <- cell_metadata.dt %>%
    .[!celltype.mapped%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}
table(cell_metadata.dt$class)
table(cell_metadata.dt$celltype.mapped)

################
## Load query ##
################

celltype_proportions_query.dt <- cell_metadata.dt %>% copy %>%
  .[,N:=.N,by="alias"] %>%
  setnames("celltype.mapped","celltype") %>%
  .[,.(celltype_proportion=.N/unique(N)),by=c("alias","class","celltype","dataset")]

sample_metadata_query.dt <- cell_metadata.dt %>%
  .[,c("alias","stage","class")] %>% unique %>%
  .[,dataset:="Query"]

################
## Load atlas ##
################

cell_metadata_atlas.dt <- fread(io$atlas.metadata) %>% 
  .[stage!="mixed_gastrulation"] %>%
  .[,sample:=as.character(sample)] %>% setnames("sample","alias")
  
sample_metadata_atlas.dt <- cell_metadata_atlas.dt %>%
  .[,c("alias","stage")] %>% unique %>%
  .[,c("class","dataset"):=list("Atlas","Atlas")]

# Load cell type proportions in the atlas
celltype_proportions_atlas.dt <- fread(paste0(io$atlas.basedir,"/results/celltype_proportions/celltype_proportions.txt.gz")) %>%
  .[sample%in%cell_metadata_atlas.dt$alias] %>%
  setnames("sample","alias") %>%
  .[,c("dataset","class"):=list("Atlas","Atlas")]


#################
## Concatenate ##
#################

sample_metadata_joint.dt <- rbind(
  sample_metadata_query.dt[,c("alias","class","stage","dataset")], 
  sample_metadata_atlas.dt[,c("alias","class","stage","dataset")]
)

celltype_proportions.dt <- rbind(
  celltype_proportions_query.dt[,c("alias","class","celltype","celltype_proportion","dataset")],
  celltype_proportions_atlas.dt[,c("alias","class","celltype","celltype_proportion","dataset")]
)

######################
## Plot atlas cells ##
######################

# celltype_proportions_atlas.mtx <- celltype_proportions_atlas.dt %>% 
#   dcast(alias~celltype, fill=0, value.var="celltype_proportion") %>% 
#   matrix.please
# 
# # PCA
# atlas.pca <- prcomp(celltype_proportions_atlas.mtx, rank.=5)
# 
# # Plot
# to.plot <- atlas.pca$x %>% as.data.table %>% 
#   .[,alias:=rownames(atlas.pca$x)] %>%
#   merge(sample_metadata_atlas.dt, by="alias")
# 
# # Flip order of PC1 such that E6.5<E8.5
# tmp <- to.plot[,mean(PC1),by="stage"]
# if (tmp[stage=="E6.5",V1]>tmp[stage=="E8.5",V1]) { to.plot[,PC1:=-PC1] }
# 
# p <- ggplot(to.plot, aes(x=PC1, y=PC2, fill=stage)) +
#   geom_point(stroke=0.5, color="black", size=4, shape=21) +
#   scale_fill_manual(values=opts$stage.colors) +
#   guides(fill=guide_legend(override.aes=list(shape=21))) +
#   guides(shape=guide_legend(override.aes=list(fill="black"))) +
#   theme_classic() +
#   theme(
#     legend.position = "right",
#     legend.title = element_blank(),
#     axis.text = element_text(color="black", size=rel(0.8))
#   )
# 
# pdf(paste0(io$outdir,"/pca_celltype_proportions_atlas.pdf"), width=7, height=5)
# print(p)
# dev.off()


##########################
## Plot cells per class ##
##########################

for (i in opts$classes) {
  
  celltype_proportions.mtx <- celltype_proportions.dt[class%in%c(i,"Atlas")] %>% 
    dcast(alias~celltype, fill=0, value.var="celltype_proportion") %>% 
    matrix.please
  
  # PCA
  pca <- prcomp(celltype_proportions.mtx, rank.=5)
  var.explained <- round(100*(pca$sdev**2)/sum(pca$sdev**2),2)
  
  # Plot
  to.plot <- pca$x %>% as.data.table %>% 
    .[,alias:=rownames(pca$x)] %>%
    merge(sample_metadata_joint.dt, by="alias")
  
  # Flip order of PC1 such that E6.5<E8.5
  tmp <- to.plot[class=="Atlas",mean(PC1),by="stage"]
  if (tmp[stage=="E6.5",V1]>tmp[stage=="E8.5",V1]) { to.plot[,PC1:=-PC1] }
  
  p1 <- ggplot(to.plot, aes(x=PC1, y=PC2, fill=stage, shape=class)) +
    geom_point(stroke=0.5, color="black", size=3.5) +
    scale_shape_manual(values=c(21,24)) +
    scale_fill_manual(values=opts$stage.colors) +
    guides(fill=guide_legend(override.aes=list(shape=21))) +
    guides(shape=guide_legend(override.aes=list(fill="black"))) +
    labs(x=sprintf("PC1 (%.2f%%)",var.explained[1]), y=sprintf("PC2 (%.2f%%)",var.explained[2])) +
    theme_classic() +
    theme(
      legend.position = "none",
      legend.title = element_blank(),
      axis.text = element_text(color="black", size=rel(0.8))
    )
  
  p2 <- ggplot(to.plot, aes(x=PC1, y=PC2, fill=class, alpha=class)) +
    geom_point(stroke=0.5, color="black", size=3.5, shape=21) +
    scale_fill_manual(values=c("gray60","red")) +
    scale_alpha_manual(values=c(0.5,1)) +
    guides(fill=guide_legend(override.aes=list(shape=21))) +
    guides(shape=guide_legend(override.aes=list(fill="black"))) +
    labs(x=sprintf("PC1 (%.2f%%)",var.explained[1]), y=sprintf("PC2 (%.2f%%)",var.explained[2])) +
    theme_classic() +
    theme(
      legend.position = "none",
      legend.title = element_blank(),
      axis.text = element_text(color="black", size=rel(0.8))
    )
  
  pdf(file.path(io$outdir,sprintf("pca_celltype_proportions_%s.pdf",i)), width=7.5, height=2.5)
  print(cowplot::plot_grid(plotlist=list(p1,p2)))
  dev.off()
}


################################
## Plot stage for each sample ##
################################

embryo_stage_list <- list()
# i <- "Dnmt1_KO"
for (i in opts$classes) {
  
  # PCA
  celltype_proportions.mtx <- celltype_proportions.dt[class%in%c(i,"Atlas")] %>% 
    dcast(alias~celltype, fill=0, value.var="celltype_proportion") %>% 
    matrix.please
  tmp <- prcomp(celltype_proportions.mtx, rank.=2)$x %>% as.data.table(keep.rownames = T) %>% 
    setnames("rn","alias") %>% merge(sample_metadata_joint.dt, by="alias")
  
  # Calculate euclidean distance between query samples and atlas samples
  foo <- tmp[class==i,c("alias","PC1","PC2")] %>% matrix.please
  bar <- tmp[class=="Atlas",c("stage","PC1","PC2")] %>% .[,.(PC1=mean(PC1), PC2=mean(PC2)), by="stage"] %>% matrix.please
  dists <- pdist::pdist(foo, bar)
  euclidean_dist.mtx <- 1/as.matrix(dists)
  rownames(euclidean_dist.mtx) <- rownames(foo); colnames(euclidean_dist.mtx) <- rownames(bar)
  euclidean_dist_norm.mtx <- sweep(euclidean_dist.mtx, 1, rowSums(euclidean_dist.mtx), "/")
  
  # Prepare data
  to.plot <- euclidean_dist_norm.mtx %>% as.data.table(keep.rownames = T) %>% setnames("rn","alias") %>%
    melt(id.vars="alias", variable.name="stage") %>% 
    .[,stage:=factor(stage,levels=opts$atlas.stages)]
  
  # Sort samples
  embryo.order <- to.plot[stage=="E8.5",mean(value),by="alias"] %>% setorder("V1") %>% .$alias
  to.plot[,alias:=factor(alias,levels=rev(embryo.order))]
  
  # Save data
  embryo_stage_list[[i]] <- to.plot %>% .[,class:=i]
  
  # Plot
  p <- ggplot(to.plot, aes(x = "", y = value, fill = stage)) + 
    geom_bar(stat = "identity", width = 1, position = position_fill(), color="black") +
    coord_polar(theta = "y") + 
    scale_fill_manual(values=opts$stage.colors) +
    # ggrepel::geom_text_repel(aes(label = stage), data=to.plot.text) +
    facet_wrap(~alias) +
    theme_classic() +
    ggplot_theme_NoAxes() +
    theme(
      legend.title = element_blank(),
      # axis.text = element_blank(),
      # axis.title = element_blank(),
      axis.line = element_blank()
      # axis.ticks = element_blank()
    )
  
  pdf(file.path(io$outdir,sprintf("pieplot_embryo_stage_%s.pdf",i)), width=11, height=8)
  print(p)
  dev.off()
  
}

###############################
## Plot stage for each class ##
###############################

embryo_stage.dt <- rbindlist(embryo_stage_list) %>% .[,class:=factor(class,levels=opts$classes)]
# fwrite(embryo_stage.dt, file.path(io$outdir,"embryo_staging.txt.gz", sep="\t", quote=F))

tmp <- embryo_stage.dt[,.(value=mean(value), sd=sd(value)), by=c("stage","class")]

p <- ggplot(embryo_stage.dt, aes_string(x="stage", y="value", fill="stage")) +
  # geom_boxplot(outlier.shape=NA, alpha=0.8) +
  geom_jitter(size=1.25, shape=21, width=0.1, alpha=0.75) +
  geom_bar(stat="identity", color="black", width=0.80, data=tmp) +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=0.25, alpha=0.75, size=0.5, data=tmp) +
  facet_wrap(~class, nrow=2) +
  # guides(x = guide_axis(angle = 90)) +
  scale_fill_manual(values=opts$stage.colors) +
  labs(x="", y="Stage probability") +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size=rel(1.20), color="black"),
    axis.text.x = element_text(size=rel(1.0), color="black"),
    axis.text.y = element_text(size=rel(1), color="black")
  )

pdf(file.path(io$outdir,"boxplots_embryo_stage_all_classes.pdf"), width=8.5, height=6.5)
print(p)
dev.off()


