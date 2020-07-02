#####################
## Define settings ##
#####################

source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
# source("/Users/ricard/10x_gastrulation_DNMTs/mapping/plot/plot_utils.R")

# io$mapping <- "/Users/ricard/data/10x_gastrulation_DNMTs/mapping/sample_metadata_mapping_mnn.txt"
io$outdir <- paste0(io$basedir,"/results/mapping/pdf")

###############
## Load data ##
###############

# mapping <- fread(io$mapping)
# sample_metadata <- sample_metadata %>% merge(mapping)

################
## Parse data ##
################

to.plot <- sample_metadata %>%
  # .[class%in%opts$classes] %>%
  .[batch%in%opts$batches] %>%
  .[!is.na(celltype.mapped),.N, by=c("stage","celltype.mapped","batch","class")] %>%
  # .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped," ", "_")] %>%
  .[,celltype.mapped:=factor(celltype.mapped,levels=opts$celltypes)]

##########
## Plot ##
##########

for (i in unique(to.plot$batch)) {
  p <- ggplot(to.plot[batch==i], aes(x=celltype.mapped, y=N)) +
    geom_bar(aes(fill=celltype.mapped), stat="identity", color="black") +
    scale_fill_manual(values=opts$celltype.colors, drop=FALSE) +
    scale_x_discrete(drop=FALSE) +
    # facet_wrap(~batch, nrow=1, scales="fixed") +
    coord_flip() +
    labs(y="Number of cells") +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=rel(1.2)),
      axis.title.x = element_text(color="black", size=rel(1.1)),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1.3), color="black"),
      axis.text.x = element_text(size=rel(1.1), color="black")
    )
  
  pdf(sprintf("%s/barplots_%s.pdf",io$outdir,i), width=9, height=7)
  print(p)
  dev.off()
}
