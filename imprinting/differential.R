##############
## Settings ##
##############

# Load settings
source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
source("/Users/ricard/10x_gastrulation_DNMTs/differential/analysis/utils.R")

# I/O
io$diff.dir <- paste0(io$basedir,"/results/differential")
io$outdir <- paste0(io$basedir,"/results/imprinting/differential")

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

# Define groups
opts$groupA <- c(
  # "E8.5_Dnmt3aHET_Dnmt3bWT",
  "E8.5_Dnmt3aKO_Dnmt3bWT",
  "E8.5_Dnmt3aKO_Dnmt3bHET",
  "E8.5_Dnmt3aWT_Dnmt3bKO",
  "E8.5_Dnmt3aHET_Dnmt3bKO",
  "E8.5_Dnmt3aKO_Dnmt3bKO",
  "E8.5_Dnmt1KO"
)

opts$groupB <- c("E8.5_WT")

opts$min.cells <- 25

#############################################
## Load results from differential analysis ##
#############################################

dt <- opts$groupA %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- sprintf("%s/%s_vs_%s_%s.txt.gz", io$diff.dir,i,opts$groupB,j)
  if (file.exists(file)) fread(file) %>% .[,c(1,2,4,6,7,10,11)] %>% .[,c("celltype","groupA","groupB"):=list(j,i,opts$groupB)]
}) %>% rbindlist }) %>% rbindlist

dt[,celltype:=factor(celltype,levels=opts$celltypes)]
dt[,groupA:=factor(groupA,levels=opts$groupA)]

# Change logFC sign (positive for upregulation events in the KOs and negative for downregulation events in the KOs)
dt[,logFC:=-logFC]

##########################
## Load imprinted genes ##
##########################

imprinting.dt <- fread("/Users/ricard/data/mm10_regulation/imprinting/parsed/mousebook_imprinted_genes.txt.gz") %>%
  setnames(c("gene","allele"))

# Subset 
dt.filt <- dt[gene%in%imprinting.dt$gene]

###############
## Bar plots ##
###############

to.plot <- dt.filt[sig==T] %>%
  .[,direction:=c("Down","Up")[as.numeric(logFC>0)+1]] %>%
  .[,.N, by=c("celltype","groupA","groupB","direction")]

p <- ggplot(to.plot, aes(x=factor(celltype), y=N)) +
  geom_bar(aes(fill = direction), color="black", stat="identity") + 
  facet_wrap(~groupA, scales="fixed") +
  labs(x="", y="Number of DE (imprinted genes)") +
  guides(x = guide_axis(angle = 90)) +
  theme_bw() +
  theme(
    legend.position = "top",
    axis.line = element_blank(),
    axis.text.x = element_text(color="black", size=rel(0.65))
  )

pdf(sprintf("%s/DE_barplots_imprinted_genes.pdf",io$outdir), width=9, height=5)
print(p)
dev.off()



###############
## Heatmaps ##
###############

to.plot <- dt.filt %>% copy# %>%
  # .[,direction:=c(-1,1)[as.numeric(logFC>0)+1]]# %>%
  # .[is.na(logFC),logFC:=0]
  # dcast(gene+groupA~celltype, value.var="logFC")

opts$max.logFC <- 3
opts$min.logFC <- -3
to.plot[logFC>opts$max.logFC,logFC:=opts$max.logFC]
to.plot[logFC<opts$min.logFC,logFC:=opts$min.logFC]

for (i in unique(to.plot$groupA)) {
  
  to.plot2 <- to.plot[groupA==i] %>%
    .[,foo:=mean(is.na(groupA_N)),by="gene"] %>% .[foo<1]
    # .[,var:=var(logFC),by="gene"] %>% .[var>0.1]
  
  p <- ggplot(to.plot2, aes(x=gene, y=celltype, fill=logFC)) +
    geom_tile() +
    # viridis::scale_fill_viridis() +
    # scale_fill_gradientn(colours = terrain.colors(10)) +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
    labs(x="", y="") +
    guides(x = guide_axis(angle = 90)) +
    theme_classic() +
    theme(
      legend.position = "top",
      axis.line = element_blank(),
      axis.text = element_text(color="black", size=rel(0.75))
    )
  
  pdf(sprintf("%s/%s_DE_heatmap_imprinted_genes.pdf",io$outdir,i), width=9, height=5)
  print(p)
  dev.off()
}

