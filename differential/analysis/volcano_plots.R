##############
## Settings ##
##############

# Load settings
source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
source("/Users/ricard/10x_gastrulation_DNMTs/differential/analysis/utils.R")

# I/O
io$diff.dir <- paste0(io$basedir,"/results/differential")
io$outdir <- paste0(io$basedir,"/results/differential/pdf/volcano_plots")

# Define groups
opts$groupA <- c(
  "E8.5_Dnmt3aKO_Dnmt3bWT",
  "E8.5_Dnmt3aHET_Dnmt3bKO",
  "E8.5_Dnmt3aHET_Dnmt3bWT",
  "E8.5_Dnmt3aKO_Dnmt3bHET",
  "E8.5_Dnmt3aKO_Dnmt3bKO",
  "E8.5_Dnmt3aWT_Dnmt3bKO",
  "E8.5_Dnmt1KO"
)
opts$groupB <- c("E8.5_WT")

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
dt <- dt[!grepl("mt-",gene)]
dt <- dt[!grepl("Rps|Rpl",gene)]
dt <- dt[!grepl("Rik",gene)]
dt <- dt[!grepl("^Hb",gene)]
dt <- dt[!gene%in%fread(io$gene_metadata)[chr=="chrY",symbol]]

# Filter by minimum number of cells per group
opts$min.cells <- 30
dt <- dt[groupA_N>opts$min.cells & groupB_N>opts$min.cells]

# Remove hits that are differentially expressed in all cell type comparisons
# foo <- dt[,mean(sig),by=c("gene")] %>% .[V1>0] %>% setorder(-V1)

##########
## Plot ##
##########

for (i in unique(dt$groupA)) {
  outdir <- sprintf("%s/%s",io$outdir,i); dir.create(outdir, showWarnings = F)
  for (j in unique(dt$celltype)) {
    
    to.plot <- dt[celltype==j & groupA==i] %>% .[!is.na(sig)] 
    if (nrow(to.plot)>0) {
      p <- gg_volcano_plot(to.plot, top_genes = 25, groupA = i, groupB = opts$groupB)
      
      pdf(sprintf("%s/%s_vs_%s_%s_volcano.pdf",outdir,i,opts$groupB,j), width=9, height=5)
      print(p)
      dev.off()
    }
  }
}
