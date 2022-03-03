source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$indir <- file.path(io$basedir,"results/differential/pseudobulk")
# io$outdir <- file.path(io$basedir,"results/differential/pseudobulk/pdf"); dir.create(io$outdir, showWarnings = F)
io$outdir <- file.path(io$basedir,"results/differential/pseudobulk/pdf"); dir.create(io$outdir, showWarnings = F)

# Options
opts$ko.classes <- c(
  "Dnmt3a_KO", 
  "Dnmt3b_KO",
  "Dnmt1_KO"
  # "Dnmt3ab_KO"
)

opts$wt.class <- "WT"

opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  # "PGC",
  # "Anterior_Primitive_Streak",
  # "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  # "Mixed_mesoderm",
  "Intermediate_mesoderm",
  "Caudal_Mesoderm",
  "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "Pharyngeal_mesoderm",
  "Cardiomyocytes",
  "Allantois",
  "ExE_mesoderm",
  "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  # "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm"
  # "ExE_ectoderm"
  # "Parietal_endoderm"
)

opts$rename_celltypes <- c(
  "Erythroid3" = "Erythroid",
  "Erythroid2" = "Erythroid",
  "Erythroid1" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors"
  # "Intermediate_mesoderm" = "Mixed_mesoderm",
  # "Paraxial_mesoderm" = "Mixed_mesoderm",
  # "Nascent_mesoderm" = "Mixed_mesoderm",
  # "Pharyngeal_mesoderm" = "Mixed_mesoderm"
  # "Visceral_endoderm" = "ExE_endoderm"
)

##############################
## Load precomputed results ##
##############################

# i <- opts$ko.classes[1]; j <- opts$celltypes[10]
diff.dt <- opts$ko.classes %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- sprintf("%s/%s/%s_%s_vs_%s.txt.gz", io$indir,i,j,opts$wt.class,i)
  if (file.exists(file)) {
    fread(file) %>% .[,c("celltype","class"):=list(j,i)]
  }
}) %>% rbindlist }) %>% rbindlist %>%
  .[,celltype:=factor(celltype,levels=opts$celltypes)] %>%
  .[,class:=factor(class,levels=opts$ko.classes)]

# Merge celltypes 
diff.dt <- diff.dt %>% 
  .[,celltype:=stringr::str_replace_all(celltype,opts$rename_celltypes)] %>%
  .[,.(expr_ko=mean(expr_ko), expr_wt=mean(expr_wt), diff=mean(diff)), by=c("celltype","class","gene")]

# save
fwrite(diff.dt, file.path(io$outdir,"diff_pseudobulk.txt.gz"))

######################
## Heatmap per gene ##
######################

# genes.to.plot <- c("Pim2", "Pou5f1", "Slc7a3", "Utf1", "Dppa5a")
# genes.to.plot <- fread(io$atlas.marker_genes)$gene %>% unique %>% .[!grepl("Rik$",.)] %>% head(n=3)
genes.to.plot <- c("Hoxc9","Hoxc8","Hoxb9","Hoxa9")

stopifnot(genes.to.plot%in%unique(diff.dt$gene))

# i <- "Hoxc8"
for (i in genes.to.plot) {
  
  to.plot <- expand.grid(X = unique(diff.dt$celltype), Y = unique(diff.dt$class)) %>% 
    as.data.table %>% setnames(c("celltype","class")) %>%
    merge(diff.dt[gene==i], by=c("celltype","class"), all.x=T) %>%
    .[,class:=factor(class,levels=rev(opts$ko.classes))]
  
  p <- ggplot(to.plot, aes(x=celltype, y=class, fill=diff)) +
    geom_tile(color="black") +
    # scale_fill_gradientn(colours = terrain.colors(10), na.value = 'gray70') +
    # scale_fill_manual(values = wes_palette("GrandBudapest1", n = 3))
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = 'gray70' ) +
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



#######################
## Heatmap all genes ##
#######################

genes.to.plot <- c("Trh","Rhox5","Rhox6","Rhox9","Trap1a","Xlr3a")
genes.to.plot <- c("Hoxc9","Hoxc8","Hoxb9","Hoxa9")
genes.to.plot <-  c("Pou5f1", "Zfp42", "Utf1", "Pim2", "Slc7a3", "Dppa5a", "Fgf5", "Dppa3", "Dppa4", "Tfap2c", "Nanog","Nr0b1","Pecam1","Gng3","Bex1","Dnmt3b")
# genes.to.plot <-  c("Bex1","Gng3","Dnmt3b")
stopifnot(genes.to.plot%in%unique(diff.dt$gene))

# i <- "Hoxc8"
to.plot <- expand.grid(X = unique(diff.dt$celltype), Y = unique(diff.dt$class)) %>% 
  as.data.table %>% setnames(c("celltype","class")) %>%
  merge(diff.dt[gene%in%genes.to.plot], by=c("celltype","class"), all.x=T) %>%
  .[,class:=factor(class,levels=opts$ko.classes)]

to.plot <- to.plot[!is.na(gene)]

p <- ggplot(to.plot, aes(x=celltype, y=class, fill=diff)) +
  geom_tile(color="black") +
  # scale_fill_gradientn(colours = terrain.colors(10), na.value = 'gray70') +
  # scale_fill_manual(values = wes_palette("GrandBudapest1", n = 3))
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = 'gray70' ) +
  facet_wrap(~gene) +
  theme_classic() +
  guides(x = guide_axis(angle = 90)) +
  theme(
    axis.text.x = element_text(color="black", size=rel(0.6)),
    axis.text.y = element_text(color="black", size=rel(0.9)),
    axis.title = element_blank(),
    strip.background = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.title = element_blank()
  )
  
  
pdf(file.path(io$outdir,"logFC_heatmap_pseudobulk_ExE_genes.pdf"), width=8, height=5)
print(p)
dev.off()
