
#######################################
## Script to plot mapping statistics ##
#######################################

library(data.table)
library(purrr)
library(ggplot2)
library(ggpubr)

source("/Users/ricard/10x_gastrulation_DNMTs/mapping/plot/plot_settings.R")

################
## Define I/O ##
################

io <- list()
# io$sample_metadata <- "/Users/ricard/data/10x_gastrulation_DNMTs/sample_metadata.txt"
io$sample_metadata <- "/Users/ricard/data/10x_gastrulation_DNMTs/mapping/sample_metadata_mapping_mnn.txt"
io$outdir <- "/Users/ricard/data/10x_gastrulation_DNMTs/mapping/pdf"

####################
## Define options ##
####################

opts <- list()

# Figure dimensions
opts$width <- 3
opts$height <- 4

###############
## Load data ##
###############

sample_metadata <- fread(io$sample_metadata) %>%
  .[,stage:="E8.5"]
  # .[!celltype.mapped.level1 %in% c("PGC","Allantois","NOIDEA",NA)]

#############################################
## Plot number of cells for each cell type ##
#############################################

foo <- sample_metadata %>%
  .[, Ncells:=.N, by=c("stage","genotype")] %>%
  .[, .(prop=.N/unique(Ncells)), by=c("celltype.mapped.level1","stage","genotype")] 

xlim.max <- max(sample_metadata[,.N,by=c("celltype.mapped.level1","stage","genotype")] %>% .[,N])


to.plot <- sample_metadata %>%
  .[,.N,by=c("celltype.mapped.level1","genotype")]

# Define colors
all(to.plot$celltype.mapped.level1 %in% names(celltype_colours))

colors.tmp <- celltype_colours[names(celltype_colours) %in% to.plot$celltype.mapped.level1]
to.plot[,celltype.mapped.level1:=factor(celltype.mapped.level1,levels=rev(names(celltype_colours)))]
to.plot[,genotype:=factor(genotype,levels=c("WT","TKO"))]

p <- barplot.pub(to.plot, x="celltype.mapped.level1", colors=celltype_colours, xlim.max=xlim.max) +
  facet_wrap(~genotype, nrow=1, scales="free_x")
  # theme(
  #   strip.background = element_blank(),
  #   strip.text = element_blank()
  # )

pdf(paste0(io$outdir,"/pdf/mapping_stats_day2.pdf"), width=opts$width, height=opts$height)
print(p)
dev.off()
