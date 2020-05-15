#####################
## Define settings ##
#####################

source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
# source("/Users/ricard/10x_gastrulation_DNMTs/mapping/plot/plot_utils.R")

io$mapping <- "/Users/ricard/data/10x_gastrulation_DNMTs/mapping/sample_metadata_mapping_mnn.txt"
io$outdir <- paste0(io$basedir,"/mapping/pdf")

###############
## Load data ##
###############

mapping <- fread(io$mapping)
sample_metadata <- sample_metadata %>% merge(mapping)

################
## Parse data ##
################

to.plot <- sample_metadata[,.N, by=c("stage","celltype.mapped.level2","genotype")] %>% 
  .[, lineage:=stringr::str_replace_all( celltype.mapped.level2,"_"," ")] %>%
  # .[, lineage:=factor( lineage,levels=names(colors))] %>%
  .[complete.cases(.)]

# QC
stopifnot(all(to.plot$lineage %in% names(colors)))

# remove ExE lineages
to.plot <- to.plot[!grepl("ExE",lineage)]

# Remove lineages with small number of cells
to.plot <- to.plot[N>100]
  
##########
## Plot ##
##########

p <- ggplot(to.plot, aes(x= lineage, y=N)) +
  geom_bar(aes(fill=lineage), stat="identity", color="black") +
  scale_fill_manual(values=colors) +
  facet_wrap(~genotype, nrow=1, scales="fixed") +
  coord_flip() +
  labs(y="Number of cells") +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(color="black", size=rel(1.3)),
    axis.title.x = element_text(color="black", size=rel(1.1)),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size=rel(1.3), color="black"),
    axis.text.x = element_text(size=rel(1.1), color="black")
  )

pdf(paste0(io$outdir,"/mapping_stats.pdf"), width=11, height=5)
print(p)
dev.off()
