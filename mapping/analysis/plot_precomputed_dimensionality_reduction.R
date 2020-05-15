###################################################################
## Plot dimensionality reduction of EB cells mapped to the atlas ##
###################################################################

# This script requires the cell metadata from the atlas, which contains the precomputed UMAP coordinates

library(data.table)
library(purrr)
library(ggplot2)

source("/Users/ricard/10x_gastrulation_DNMTs/mapping/plot/plot_utils.R")

#########
## I/O ##
#########

io <- list()
io$path2atlas <- "/Users/ricard/data/gastrulation10x"
io$path2query <- "/Users/ricard/data/10x_gastrulation_DNMTs"
io$mapping <- "/Users/ricard/data/10x_gastrulation_DNMTs/mapping/mapping10x_mnn.rds"
io$outdir <- "/Users/ricard/data/10x_gastrulation_DNMTs/mapping/pdf"

#############
## Options ##
#############

opts <- list()

# Dot size
opts$size.mapped <- 0.5
opts$size.nomapped <- 0.1

# Transparency
opts$alpha.mapped <- 1.0
opts$alpha.nomapped <- 0.35

####################
## Load 10x atlas ##
####################

# Load atlas cell metadata
meta_atlas <- fread(paste0(io$path2atlas, "/sample_metadata.txt")) %>%
  .[stripped==F & doublet==F]

# Extract precomputed dimensionality reduction coordinates
umap <- meta_atlas[,c("cell","umapX","umapY","celltype")] %>%
  setnames(c("umapX","umapY"),c("V1","V2"))

#####################
## Load query data ##
#####################

# Load query cell metadata
meta_query <- fread(paste0(io$path2query, "/mapping/sample_metadata_mapping_mnn.txt"))

# Load precomputed mapping
mapping <- readRDS(io$mapping)$mapping
mapping.dt <- data.table(
  cell            = mapping$cell, 
  celltype.mapped = mapping$celltype.mapped,
  stage.mapped    = mapping$stage.mapped,
  closest.cell    = as.character(mapping$closest.cell)
)

################
## Parse data ##
################

# Remove lineages
mapping.dt <- mapping.dt[!celltype.mapped%in%c("ExE ectoderm","PGC", "Parietal endoderm")]
umap <- umap[!celltype%in%c("ExE ectoderm","PGC", "Parietal endoderm")]

###################################
## Plot dimensionality reduction ##
###################################

# Prepare query data.frame to plot
plot_df_query = mapping.dt %>% merge(meta_query, by="cell")

# Prepare atlas data.frame to plot
plot_df_atlas = umap %>% merge(meta_atlas, by="cell")

plot_df_atlas[,index.wt:=match(plot_df_atlas$cell, plot_df_query[genotype=="WT",closest.cell] )]
plot_df_atlas[,index.ko:=match(plot_df_atlas$cell, plot_df_query[genotype=="TKO",closest.cell] )]
plot_df_atlas[,mapped.wt:=c(0,-10)[as.numeric(as.factor(!is.na(index.wt)))]]
plot_df_atlas[,mapped.ko:=c(0,10)[as.numeric(as.factor(!is.na(index.ko)))]]
plot_df_atlas[,mapped:=factor(mapped.wt + mapped.ko, levels=c("0","-10","10"))] %>% setorder(mapped)

p <- plot.dimred.wtko(plot_df_atlas) +
  theme(legend.position = "top", legend.title = element_blank())

pdf(paste0(io$outdir,"/umap_mapped.pdf"), width=8, height=6.5)
print(p)
dev.off()
