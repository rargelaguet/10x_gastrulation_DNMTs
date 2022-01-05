source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$indir <- file.path(io$basedir,"results_new/differential/pseudobulk")
io$outdir <- file.path(io$basedir,"results_new/differential/pseudobulk/pdf/individual_genes"); dir.create(io$outdir, showWarnings = F)

# Options
opts$ko.classes <- c(
  "E8.5_Dnmt3aKO_Dnmt3bWT",
  # "E8.5_Dnmt3aHET_Dnmt3bKO",
  # "E8.5_Dnmt3aHET_Dnmt3bWT",
  # "E8.5_Dnmt3aKO_Dnmt3bHET",
  "E8.5_Dnmt3aKO_Dnmt3bKO",
  "E8.5_Dnmt3aWT_Dnmt3bKO",
  "E8.5_Dnmt1KO"
)

opts$wt.class <- "E8.5_WT"

##############################
## Load precomputed results ##
##############################

i <- opts$ko.classes[1]
diff.dt <- opts$ko.classes %>% map(function(i) {
  file <- sprintf("%s/%s/%s_vs_%s.txt.gz", io$indir,i,opts$wt.class,i)
  if (file.exists(file)) {
    fread(file) %>% .[,c("class"):=list(i)]
  }
}) %>% rbindlist %>%
  .[,class:=factor(class,levels=opts$ko.classes)]

######################
## Heatmap per gene ##
######################

# genes.to.plot <- c("Pim2", "Pou5f1", "Slc7a3", "Utf1", "Dppa5a")
genes.to.plot <- fread(io$atlas.marker_genes)$gene %>% unique %>% .[!grepl("Rik$",.)]

ylimits <- 4

# i <- "Hoxc8"
for (i in genes.to.plot) {
  
  to.plot <- diff.dt[gene==i]
  
  p <- ggplot(to.plot, aes(x=class, y=diff)) +
    # geom_bar(stat="identity", color="black") +
    # geom_segment(xend=0, yend=0) +
    geom_point(shape=21) +
    coord_flip(ylim=c(-ylimits,ylimits)) +
    geom_hline(yintercept=0, linetype="dashed", size=0.5) +
    # scale_fill_gradientn(colours = terrain.colors(10), na.value = 'gray70') +
    # scale_fill_manual(values = wes_palette("GrandBudapest1", n = 3))
    # scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = 'gray70' ) +
    theme_classic() +
    guides(x = guide_axis(angle = 90)) +
    theme(
      axis.text = element_text(color="black", size=rel(0.7)),
      axis.title = element_blank(),
      strip.background = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      legend.title = element_blank()
    )
  
  
  pdf(file.path(io$outdir,sprintf("%s_logFC_heatmap_pseudobulk.pdf",i)), width=8, height=5)
  print(p)
  dev.off()
}
