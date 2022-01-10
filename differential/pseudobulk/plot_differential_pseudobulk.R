source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$indir <- file.path(io$basedir,"results_all/differential/pseudobulk")
io$outdir <- file.path(io$basedir,"results_all/differential/pseudobulk/pdf/individual_genes"); dir.create(io$outdir, showWarnings = F)

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

opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  # "PGC",
  "Anterior_Primitive_Streak",
  "Notochord",
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
genes.to.plot <- fread(io$atlas.marker_genes)$gene %>% unique %>% .[!grepl("Rik$",.)]

# i <- "Hoxc8"
for (i in genes.to.plot) {
  
  to.plot <- expand.grid(X = unique(diff.dt$celltype), Y = unique(diff.dt$class)) %>% 
    as.data.table %>% setnames(c("celltype","class")) %>%
    merge(diff.dt[gene==i], by=c("celltype","class"), all.x=T)
  
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
