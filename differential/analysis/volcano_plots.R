##############
## Settings ##
##############

# Load settings
source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
source("/Users/ricard/10x_gastrulation_DNMTs/differential/analysis/utils.R")

# I/O
io$diff.dir <- paste0(io$basedir,"/results/differential")
io$outdir <- paste0(io$basedir,"/results/differential/pdf")

# Define groups
opts$groupA <- c(
  # "E8.5_Dnmt3aKO_Dnmt3bWT"
  # "E8.5_Dnmt3aHET_Dnmt3bKO", 
  # "E8.5_Dnmt3aHET_Dnmt3bWT", 
  # "E8.5_Dnmt3aKO_Dnmt3bHET", 
  "E8.5_Dnmt3aKO_Dnmt3bKO"
  # "E8.5_Dnmt3aWT_Dnmt3bKO"
)
opts$groupB <- c("E8.5_Dnmt3aWT_Dnmt3bWT")

opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  # "Caudal_epiblast",
  # "PGC",
  # "Anterior_Primitive_Streak",
  # "Notochord",
  # "Def._endoderm",
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

#############################################
## Load results from differential analysis ##
#############################################

dt <- opts$groupA %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- sprintf("%s/%s_vs_%s_%s.txt.gz", io$diff.dir,i,opts$groupB,j)
  if (file.exists(file)) fread(file) %>% .[,c(1,2,4,6,7,10,11)] %>% 
    setnames(c(sprintf("N_%s",i),sprintf("N_%s",opts$groupB)),c("N_groupA","N_groupB")) %>%
    .[,c("celltype","groupA","groupB"):=list(j,i,opts$groupB)]
}) %>% rbindlist }) %>% rbindlist

####################
## Filter results ##
####################

# Remove some hits
dt <- dt[gene!="Xist"]

# Filter by minimum number of cells per group
opts$min.cells <- 30
dt <- dt[N_groupA>opts$min.cells & N_groupB>opts$min.cells]

# Remove manual hits
dt <- dt[gene!="Xist"]

# Remove hits that are differentially expressed in all cell type comparisons
# foo <- dt[,mean(sig),by=c("gene")] %>% .[V1>0] %>% setorder(-V1)

##########
## Plot ##
##########

for (i in unique(dt$celltype)) {
  to.plot <- dt[celltype==i] %>% .[!is.na(sig)] 
  p <- gg_volcano_plot(to.plot, top_genes=20)
  
  pdf(sprintf("%s/%s_vs_%s_%s_volcano.pdf",io$outdir,opts$groupA,opts$groupB,i), width=9, height=5)
  print(p)
  dev.off()
}
