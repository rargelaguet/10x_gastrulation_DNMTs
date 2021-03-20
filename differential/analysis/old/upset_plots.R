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
# "E8.5_Dnmt3aWT_Dnmt3bWT",

opts$groupB <- c("E8.5_Dnmt3aWT_Dnmt3bWT")

#############################################
## Load results from differential analysis ##
#############################################

dt <- opts$groupA %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- sprintf("%s/%s_vs_%s_%s.txt.gz", io$diff.dir,i,opts$groupB,j)
  if (file.exists(file)) fread(file) %>% .[,c(1,2,4,10,11)] %>% .[,c("celltype","groupA","groupB"):=list(j,i,opts$groupB)]
}) %>% rbindlist }) %>% rbindlist

unique(dt$groupA)
unique(dt$groupB)
unique(dt$celltype)

# Remove some hits
dt <- dt[gene!="Xist"]

###########
## Parse ##
###########

dt[is.na(sig),sig:=F]
sum(is.na(dt$sig))

foo <- dt[,mean(sig),by=c("gene","groupA")] %>% 
  setorder(-V1) %>% .[V1>0]

foo[groupA=="E8.5_Dnmt3aKO_Dnmt3bKO"] %>% View

dt[sig==T & grepl("Hox",gene) & groupA=="E8.5_Dnmt3aKO_Dnmt3bKO"] %>% View

##########
## Plot ##
##########

