
here::i_am("gene_regulatory_networks/gene_regulatory_networks.R")


source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O ##
io$grn <- "/Users/argelagr/data/gastrulation_multiome_10x/results_new/rna_atac/gene_regulatory_networks/pseudobulk/TF2gene_after_virtual_chip.txt.gz"
io$outdir <- file.path(io$basedir,"results_new/individual_genes")

## Define options ##

# Define cell types to plot
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
	"Blood_progenitors",
	# "Blood_progenitors_1",
	# "Blood_progenitors_2",
	"Erythroid",
	# "Erythroid1",
	# "Erythroid2",
	# "Erythroid3",
	"NMP",
	"Rostral_neurectoderm",
	"Caudal_neurectoderm",
	"Neural_crest",
	"Forebrain_Midbrain_Hindbrain",
	"Spinal_cord",
	"Surface_ectoderm"
	# "Visceral_endoderm",
	# "ExE_endoderm",
	# "ExE_ectoderm",
	# "Parietal_endoderm"
)

# Define classes to plot
opts$classes <- c(
  "E8.5_WT"
  # "E8.5_Dnmt3aHET_Dnmt3bWT",
  # "E8.5_Dnmt3aKO_Dnmt3bWT",
  # "E8.5_Dnmt3aKO_Dnmt3bHET",
  # "E8.5_Dnmt3aWT_Dnmt3bKO",
  # "E8.5_Dnmt3aHET_Dnmt3bKO",
  # "E8.5_Dnmt3aKO_Dnmt3bKO",
  # "E8.5_Dnmt1KO"
)

opts$rename_celltypes <- c(
  "Erythroid3" = "Erythroid",
  "Erythroid2" = "Erythroid",
  "Erythroid1" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors",
  "Anterior_Primitive_Streak" = "Primitive_Streak"
)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,opts$rename_celltypes)] %>%
  .[pass_rnaQC==TRUE & celltype.mapped%in%opts$celltypes & class%in%opts$classes]

# Only consider cell types with sufficient observations in WT cells
celltypes.to.use <- sample_metadata[class=="E8.5_WT",.(N=.N),by="celltype.mapped"] %>% .[N>=50,celltype.mapped]
sample_metadata <- sample_metadata[celltype.mapped%in%celltypes.to.use]

table(sample_metadata$class)
table(sample_metadata$celltype.mapped)

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- load_SingleCellExperiment(io$sce, cells=sample_metadata$cell, normalise = TRUE)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

##############
## Load GRN ##
##############

opts$min.chip.score <- 0.25

grn.dt <- fread(io$grn) %>%
  .[,.SD[which.min(dist)],by=c("tf","gene")] %>%
  .[chip_score>=opts$min.chip.score,c("tf","chip_score","gene","dist")]

#################
## Filter data ##
#################

genes <- intersect(rownames(sce),unique(grn.dt$gene))
sce.genes <- sce[genes,]
grn.dt <- grn.dt[gene%in%genes]

TFs <- intersect(toupper(rownames(sce)), unique(grn.dt$tf))
sce.tf <- sce[stringr::str_to_title(TFs)]; rownames(sce.tf) <- toupper(rownames(sce.tf))
grn.dt <- grn.dt[tf%in%TFs]

# grn.mtx <- grn.dt %>% dcast(tf~gene,value.var="chip_score") %>% matrix.please
# grn.mtx[!is.na(grn.mtx)] <- 1

grn.list <- grn.dt %>% split(.$tf) %>% map(~ .[["gene"]])

#####################
## Mean expression ##
#####################

i <- "TAL1"
for (i in TFs) {
  sce.tmp <- sce.genes[grn.dt[tf==i,gene],]
  
  dt <- data.table(
    tf_module = i,
    expr = colMeans(logcounts(sce.tmp)),
    cell = colnames(sce.tmp)
  ) %>% merge(sample_metadata[,c("cell","class","celltype.mapped")], by="cell")
  
  ggboxplot(dt, x="celltype.mapped", y="expr", fill="celltype.mapped") + 
    scale_fill_manual(values=opts$celltype.colors) +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
}

#####################################
## Differential expression results ##
#####################################

io$indir <- file.path(io$basedir,"results_new/differential")

source(here::here("differential/analysis/load_data.R"))

i <- "TAL1"
for (i in TFs) {
  
  diff_tmp.dt <- diff.dt[class=="E8.5_Dnmt1KO" & gene%in%grn.list[["TAL1"]]]
  
  to.plot <- diff_tmp.dt[,sum(sig,na.rm=T),by="celltype"] %>% .[,tf_module:=i]
  
  ggbarplot(to.plot, x="celltype", y="V1", fill="celltype") + 
    scale_fill_manual(values=opts$celltype.colors) +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
}

############
## AUCell ##
############

# library(AUCell)
# logcounts(sce.t)
# cells_rankings <- AUCell_buildRankings(logcounts(sce.genes), nCores=1, plotStats=TRUE)
# cells_AUC <- AUCell_calcAUC(grn.list, cells_rankings)
# 
# min(getAUC(cells_AUC))
# max(getAUC(cells_AUC))
# 
# i <- "TAL1"
# for (i in TFs) {
#   
#   sce.tmp <- sce.genes[grn.dt[tf==i,gene],]
#   
#   to.plot <- data.table(
#     tf_module = i,
#     auc = getAUC(cells_AUC)[i,],
#     cell = colnames(cells_AUC)
#   ) %>% merge(sample_metadata[,c("cell","class","celltype.mapped")], by="cell")
#   
#   ggboxplot(to.plot, x="celltype.mapped", y="auc", fill="celltype.mapped", outlier.shape=NA) + 
#     scale_fill_manual(values=opts$celltype.colors) +
#     theme(
#       legend.position = "none",
#       axis.text.x = element_blank(),
#       axis.ticks.x = element_blank()
#     )
# }

##########
## Plot ##
##########

genes.to.plot <- c("Tex19.1","Morc1","Dppa3","Rex1","Dppa5a","Dppa4","Dppa2","Zfp981")
# genes.to.plot <- c("Lefty1","Cd34","Tmsb4x","Fgf3","Spata7","Cer1","Spink1","Dppa4","Dppa5a","Prc1","Lefty2","Ube2c","Hba-x","Hbb-y","Hba-a1","Hbb-bh1")
# genes.to.plot <- c("Vegfa","Vegfb","Vegfc","Vegfd","Kdr","Flt1","Tal1","Runx1","Etv2)
# genes.to.plot <- c("Tet1","Tet2","Tet3","Dnmt1","Dnmt3a","Dnmt3b","Dnmt3l")
genes.to.plot <- rownames(sce)[grep("Tet",rownames(sce))]
# genes.to.plot <- fread(io$atlas.marker_genes)$gene %>% unique %>% .[!grepl("Rik$",.)]
# genes.to.plot <- fread("/Users/ricard/data/gastrulation10x/results/differential/celltypes/E8.5/Neural_crest_vs_Forebrain_Midbrain_Hindbrain.txt.gz") %>% .[sig==T & logFC<(-2),gene]

for (i in 1:length(genes.to.plot)) {
  
  gene <- genes.to.plot[i]
  
  if (gene %in% rownames(sce)) {
    print(sprintf("%s/%s: %s",i,length(genes.to.plot),gene))
    outfile <- sprintf("%s/%s.pdf",io$outdir,gene)
  
    to.plot <- data.table(
      cell = colnames(sce),
      expr = logcounts(sce)[gene,]
    ) %>% merge(sample_metadata[,c("cell","sample","class","celltype.mapped")], by="cell") %>%
      .[,N:=.N,by=c("sample","celltype.mapped")] %>% .[N>=10]
    
    p <- ggplot(to.plot, aes(x=class, y=expr, fill=class)) +
      geom_violin(scale = "width", alpha=0.8) +
      geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.8) +
      stat_summary(fun.data = give.n, geom = "text", size=3) +
      # geom_jitter(size=2, shape=21, stroke=0.2, alpha=0.5) +
      # scale_fill_manual(values=opts$colors) +
      scale_fill_brewer(palette="Dark2") +
      facet_wrap(~celltype.mapped, scales="fixed") +
      theme_classic() +
      labs(x="",y=sprintf("%s expression",gene)) +
      guides(x = guide_axis(angle = 90)) +
      theme(
        strip.text = element_text(size=rel(0.85)),
        axis.text.x = element_text(colour="black",size=rel(0.9)),
        axis.text.y = element_text(colour="black",size=rel(0.9)),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(colour="black",size=rel(1.0)),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=rel(0.85))
      )
    
      pdf(outfile, width=10, height=9)
      print(p)
      dev.off()
        
  } else {
    print(sprintf("%s not found",gene))
  }
}




