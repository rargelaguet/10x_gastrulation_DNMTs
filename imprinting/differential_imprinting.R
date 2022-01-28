source(here::here("settings.R"))
source(here::here("utils.R"))


##############
## Settings ##
##############

# Load settings
# source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
source(here::here("differential/analysis/utils.R"))

# I/O
io$imprinted.genes <- "/Users/argelagr/data/mm10_regulation/imprinting/parsed/mousebook_imprinted_genes.txt.gz"
io$indir <- paste0(io$basedir,"/results_all/differential")
io$outdir <- paste0(io$basedir,"/results_all/imprinting")

# Options
# opts$classes <- c(
#   # "Dnmt3a_HET_Dnmt3b_KO",
#   # "Dnmt3a_HET_Dnmt3b_WT",
#   # "Dnmt3a_KO_Dnmt3b_HET",
#   "WT",
#   "Dnmt3a_KO",
#   "Dnmt3b_KO",
#   "Dnmt1_KO",
#   "Dnmt3ab_KO"
# )
opts$min.cells <- 30

#############################################
## Load results from differential analysis ##
#############################################

source(here::here("differential/analysis/load_data.R"))

##########################
## Load imprinted genes ##
##########################

imprinting.dt <- fread(io$imprinted.genes) %>%
  setnames(c("gene","allele"))

####################
## Filter results ##
####################

# Subset imprinted genes
diff_imprinted.dt <- diff.dt[gene%in%imprinting.dt$gene]

diff_imprinted.dt <- diff_imprinted.dt[groupA_N>=opts$min.cells & groupB_N>=opts$min.cells]

###############
## Bar plots ##
###############

to.plot <- diff_imprinted.dt[sig==T] %>%
  .[,direction:=c("Down","Up")[as.numeric(logFC>0)+1]] %>%
  .[,.N, by=c("celltype","class","direction")]

p <- ggplot(to.plot, aes(x=celltype, y=N)) +
  geom_bar(aes(fill = direction), color="black", stat="identity", position="dodge") + 
  facet_wrap(~class, scales="fixed") +
  labs(x="", y="Number of DE (imprinted genes)") +
  guides(x = guide_axis(angle = 90)) +
  theme_bw() +
  theme(
    legend.position = "top",
    axis.line = element_blank(),
    axis.text.x = element_text(color="black", size=rel(0.65))
  )

pdf(sprintf("%s/DE_barplots_imprinted_genes.pdf",io$outdir), width=9, height=5)
print(p)
dev.off()



###################################
## Heatmaps, one class at a time ##
###################################

opts$max.logFC <- 4; opts$min.logFC <- -4

for (i in opts$ko.classes) {
  
  to.plot <- diff_imprinted.dt[class==i] %>%
    .[logFC>opts$max.logFC,logFC:=opts$max.logFC] %>%
    .[logFC<opts$min.logFC,logFC:=opts$min.logFC]
  
  p <- ggplot(to.plot, aes(x=gene, y=celltype, fill=logFC)) +
    geom_tile() +
    # viridis::scale_fill_viridis() +
    # scale_fill_gradientn(colours = terrain.colors(10)) +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, limits=c(opts$min.logFC-0.01,opts$max.logFC+0.01)) +
    labs(x="", y="") +
    guides(x = guide_axis(angle = 90)) +
    theme_classic() +
    theme(
      legend.position = "top",
      axis.line = element_blank(),
      axis.text = element_text(color="black", size=rel(0.75))
    )
  
  pdf(sprintf("%s/%s_DE_heatmap_imprinted_genes.pdf",io$outdir,i), width=9, height=5)
  print(p)
  dev.off()
}

###########################
## Heatmaps, all classes ##
###########################

stop("FILL IN MISSING VALUES")

to.plot <- diff_imprinted.dt %>% copy %>%
  .[logFC>opts$max.logFC,logFC:=opts$max.logFC] %>%
  .[logFC<opts$min.logFC,logFC:=opts$min.logFC]

p <- ggplot(to.plot, aes(x=gene, y=celltype, fill=logFC)) +
  geom_tile(color="black") +
  facet_wrap(~class, scales="free") +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, limits=c(opts$min.logFC-0.01,opts$max.logFC+0.01)) +
  labs(x="", y="") +
  guides(x = guide_axis(angle = 90)) +
  theme_classic() +
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    axis.line = element_blank(),
    axis.text.y = element_text(color="black", size=rel(0.45)),
    axis.text.x = element_text(color="black", size=rel(0.65)),
    strip.text = element_text(color="black", size=rel(1)),
    axis.ticks = element_line(size=rel(0.75))
  )

pdf(file.path(io$outdir,"DE_heatmap_imprinted_genes.pdf"), width=15, height=9)
print(p)
dev.off()
