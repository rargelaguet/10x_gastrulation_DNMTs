##############
## Settings ##
##############

# Load settings
source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")

# I/O
io$diff.dir <- paste0(io$basedir,"/results/differential")
io$outdir <- paste0(io$basedir,"/results/differential/cellfate_bias")

opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  "PGC",
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
  "Visceral_endoderm"
  # "ExE_endoderm",
  # "ExE_ectoderm",
  # "Parietal_endoderm"
)

opts$ExE.celltypes <- c(
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm",
  "Parietal_endoderm"
)

# Define groups
opts$groupA <- c(
  # "E8.5_Dnmt3aKO_Dnmt3bWT",
  # "E8.5_Dnmt3aHET_Dnmt3bKO",
  # "E8.5_Dnmt3aHET_Dnmt3bWT",
  # "E8.5_Dnmt3aWT_Dnmt3bKO",
  # "E8.5_Dnmt3aKO_Dnmt3bHET",
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

unique(dt$groupA)
unique(dt$groupB)
unique(dt$celltype)

# dt[celltype=="Neural_crest" & sig==TRUE] %>% .[,N:=.N,by="gene"] %>% .[N>1] %>% View
# dt.filt[groupA=="E8.5_Dnmt3aKO_Dnmt3bWT" & celltype=="Caudal_epiblast" & sig==TRUE] %>% View

####################
## Filter results ##
####################

# Filter
dt.filt <- dt[groupA_N>=opts$min.cells & groupB_N>=opts$min.cells]

# Remove sex-specific differences
# dt.filt <- dt[gene!="Xist"]

# Remove Riks
dt.filt <- dt.filt[!grep("Rik",gene)]

# Remove hemoglobins
dt.filt <- dt.filt[!grep("^Hb",gene)]

# Subset to lineage markers
marker_genes.dt <- fread(io$atlas.marker_genes) %>% .[,ExE:=celltype%in%opts$ExE.celltypes]
dt.filt <- dt.filt[gene%in%unique(marker_genes.dt$gene)]
marker_genes.dt[,.N,by="celltype"]

# dt.filt[celltype=="Blood_progenitors_1" & groupA=="E8.5_Dnmt1KO"] %>%
#   .[gene%in%marker_genes.dt[celltype=="Rostral_neurectoderm",gene]]

##################################
## Quantify fate disruption (?) ##
##################################

to.plot <- dt.filt %>%
  merge(
    marker_genes.dt[,c("gene","celltype","ExE")] %>% setnames("celltype","celltype_marker"),
    by = c("gene"), allow.cartesian=TRUE
  ) %>%
  .[,mean(sig), by=c("celltype","celltype_marker","groupA","groupB")]

p <- ggplot(to.plot, aes(x=factor(celltype), y=N)) +
  geom_bar(aes(fill = celltype_marker), color="black", stat="identity") + 
  facet_wrap(~groupA, scales="free_y") +
  scale_fill_manual(values=opts$celltype.colors, drop=F) +
  labs(x="", y="Number of DE genes") +
  guides(x = guide_axis(angle = 90)) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.line = element_blank(),
    axis.text.x = element_text(color="black", size=rel(0.65))
  )

pdf(sprintf("%s/DE_barplots_marker_genes_fate_bias.pdf",io$outdir), width=9, height=5)
print(p)
dev.off()
