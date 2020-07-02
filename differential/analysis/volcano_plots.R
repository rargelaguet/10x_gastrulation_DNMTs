##############
## Settings ##
##############

source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
source("/Users/ricard/10x_gastrulation_DNMTs/differential/analysis/utils.R")
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

#############################################
## Load results from differential analysis ##
#############################################

dt <- opts$groupA %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- sprintf("%s/%s_vs_%s_%s.txt.gz", io$diff.dir,i,opts$groupB,j)
  if (file.exists(file)) fread(file) %>% .[,c(1,2,4,6,7,10,11)] %>% 
    setnames(c(sprintf("N_%s",i),sprintf("N_%s",opts$groupB)),c("N_groupA","N_groupB")) %>%
    .[,c("celltype","groupA","groupB"):=list(j,i,opts$groupB)]
}) %>% rbindlist }) %>% rbindlist

# Remove some hits
dt <- dt[gene!="Xist"]

# Filter by minimum number of cells per group
opts$min.cells <- 30
dt <- dt[N_groupA>opts$min.cells & N_groupB>opts$min.cells]

# Remove manual hits
dt <- dt[gene!="Xist"]

# Remove hits that are differentially expressed in all cell type comparisons
# foo <- dt[,mean(sig),by=c("gene")]

for (i in unique(dt$celltype)) {
  to.plot <- dt[celltype==i] %>% .[!is.na(sig)] 
  p <- gg_volcano_plot(to.plot, top_genes=20)
  
  pdf(sprintf("%s/%s_vs_%s_%s_volcano.pdf",io$outdir,opts$groupA,opts$groupB,i), width=9, height=5, useDingbats = F)
  print(p)
  dev.off()
}
