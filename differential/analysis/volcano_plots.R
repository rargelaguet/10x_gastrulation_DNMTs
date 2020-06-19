##############
## Settings ##
##############

source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
source("/Users/ricard/10x_gastrulation_DNMTs/differential/analysis/utils.R")
io$diff.dir <- paste0(io$basedir,"/results/differential")
io$outdir <- paste0(io$basedir,"/results/differential/pdf")

# Define groups
opts$groupA <- c(
  "E8.5_Dnmt3aKO_Dnmt3bWT", 
  "E8.5_Dnmt3aHET_Dnmt3bKO", 
  "E8.5_Dnmt3aHET_Dnmt3bWT", 
  "E8.5_Dnmt3aKO_Dnmt3bHET", 
  "E8.5_Dnmt3aKO_Dnmt3bKO", 
  "E8.5_Dnmt3aWT_Dnmt3bKO"
)
opts$groupB <- c("E8.5_Dnmt3aWT_Dnmt3bWT")

#############################################
## Load results from differential analysis ##
#############################################

dt <- opts$celltypes %>% map(function(i) {
  file <- sprintf("%s/%s_vs_%s_%s.txt.gz", io$diff.dir,opts$groupA,opts$groupB,i)
  if (file.exists(file)) fread(file) %>% .[,c("celltype","groupA","groupB"):=list(i,opts$groupA,opts$groupB)]
}) %>% rbindlist
# (TO-DO) Filter by minimum number of cells per group

for (i in unique(dt$celltype)) {
  to.plot <- dt[celltype==i] %>% .[!is.na(sig)] 
  p <- gg_volcano_plot(to.plot, top_genes=15)
  
  pdf(sprintf("%s/%s_vs_%s_%s_volcano.pdf",io$outdir,opts$groupA,opts$groupB,i), width=9, height=5, useDingbats = F)
  print(p)
  dev.off()
}
