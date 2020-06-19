###################################################################
## Plot dimensionality reduction of EB cells mapped to the atlas ##
###################################################################

source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
source("/Users/ricard/10x_gastrulation_DNMTs/mapping/analysis/plot_utils.R")

io$outdir <- paste0(io$basedir,"/results/mapping/pdf")

#############
## Options ##
#############

# opts$batches <- c(
#   "E8.5_Dnmt3aKO_Dnmt3bWT", 
#   "E8.5_Dnmt3aHET_Dnmt3bKO", 
#   "E8.5_Dnmt3aHET_Dnmt3bWT", 
#   "E8.5_Dnmt3aKO_Dnmt3bHET", 
#   "E8.5_Dnmt3aKO_Dnmt3bKO", 
#   "E8.5_Dnmt3aWT_Dnmt3bKO"
# )
# # "E8.5_Dnmt3aWT_Dnmt3bWT",

opts$batches <- c(
  "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001", 
  # "SIGAB6_E85_3_Dnmt3aWT_Dnmt3b_WT_L002", 
  "SIGAC6_E85_5_Dnmt3aKO_Dnmt3b_Het_L003", 
  "SIGAD6_E85_8_Dnmt3aHet_Dnmt3b_KO_L004",
  # "15_E8_5_D3A_WT_D3B_WT_L007",
  "17_E8_5_D3A_KO_D3B_WT_L008",
  "2_E8_5_D3A_WT_D3B_KO_L003",
  "3_E8_5_D3A_HET_D3B_WT_L004",
  "7_E8_5_D3A_WT_D3B_KO_L005",
  "8_E8_5_D3A_KO_D3B_KO_L006"
)

# Dot size
opts$size.mapped <- 0.3
opts$size.nomapped <- 0.1

# Transparency
opts$alpha.mapped <- 0.9
opts$alpha.nomapped <- 0.35

####################
## Load 10x atlas ##
####################

# Load atlas cell metadata
meta_atlas <- fread(io$atlas.metadata) %>%
  .[stripped==F & doublet==F]

# Extract precomputed dimensionality reduction coordinates
umap.dt <- meta_atlas[,c("cell","umapX","umapY","celltype")] %>%
  setnames(c("umapX","umapY"),c("V1","V2"))

########################################################
## Plot dimensionality reduction: one batch at a time ##
########################################################

for (i in opts$batches) {
  to.plot <- umap.dt %>% copy %>%
    .[,index:=match(cell, sample_metadata[batch==i,closest.cell] )] %>% 
    .[,mapped:=as.factor(!is.na(index))] %>% 
    .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas",i))] %>%
    setorder(mapped) 
  
  p <- plot.dimred(to.plot, query.label = i, atlas.label = "Atlas")
  
  pdf(sprintf("%s/umap_mapped_%s.pdf",io$outdir,i), width=8, height=6.5)
  print(p)
  dev.off()
}

######################################################
## Plot dimensionality reduction: WT vs KO together ##
######################################################

for (i in opts$batches) {
  to.plot <- umap.dt %>% copy %>%
    .[,index.wt:=match(cell, sample_metadata[class=="E8.5_Dnmt3aWT_Dnmt3bWT",closest.cell] )] %>%
    .[,index.ko:=match(cell, sample_metadata[batch==i,closest.cell] )] %>%
    .[,mapped.wt:=c(0,-10)[as.numeric(as.factor(!is.na(index.wt)))]] %>%
    .[,mapped.ko:=c(0,10)[as.numeric(as.factor(!is.na(index.ko)))]] %>%
    .[,mapped:=factor(mapped.wt + mapped.ko, levels=c("0","-10","10"))] %>%
    .[,mapped:=plyr::mapvalues(mapped, from = c("0","-10","10"), to = c("Atlas","WT",i))] %>% setorder(mapped)
  
  p <- plot.dimred.wtko(to.plot, wt.label = "WT", ko.label = i, nomapped.label = "Atlas")
  
  pdf(sprintf("%s/umap_mapped_%s.pdf",io$outdir,i), width=8, height=6.5)
  print(p)
  dev.off()
}
