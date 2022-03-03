source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$indir <- file.path(io$basedir,"results/differential/pseudobulk")
io$outdir <- file.path(io$basedir,"results/differential/pseudobulk/pdf/fig"); dir.create(io$outdir, showWarnings = F, recursive = T)

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
)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  setnames("celltype.mapped","celltype") %>%
  .[pass_rnaQC==TRUE & celltype%in%opts$celltypes & class%in%c("WT",opts$ko.classes)] %>%
  .[,celltype:=stringr::str_replace_all(celltype,opts$rename_celltypes)]

celltype_numbers.dt <- sample_metadata[,.(Ncells=.N, Nembryos=length(unique(sample))), by=c("class","celltype")]

# Require at least 50 cells and 3 samples per condition
celltype_numbers.dt <- celltype_numbers.dt[Nembryos>=3 & Ncells>=50]

##############################
## Load precomputed results ##
##############################

diff.dt <- opts$ko.classes %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- sprintf("%s/%s/%s_%s_vs_%s.txt.gz", io$indir,i,j,opts$wt.class,i)
  if (file.exists(file)) {
    fread(file) %>% .[,c("celltype","class"):=list(j,i)]
  }
}) %>% rbindlist }) %>% rbindlist

# Merge celltypes 
diff.dt <- diff.dt %>% 
  .[,celltype:=stringr::str_replace_all(celltype,opts$rename_celltypes)] %>%
  .[,.(expr_ko=mean(expr_ko), expr_wt=mean(expr_wt), diff=round(mean(diff),2)), by=c("celltype","class","gene")]

#######################
## Heatmap ExE genes ##
#######################

genes.to.plot <- c("Rhox9","Trap1a","Xlr3a","Ascl2","Apoe","Fmr1nb","Tex19.1") # "Trh","Rhox5","Rhox6"
# celltypes.to.plot <- c("Caudal_epiblast", "Gut", "Cardiomyocytes", "Rostral_neurectoderm", "Paraxial_mesoderm", "Spinal_cord")
# celltypes.to.plot <- c("Cardiomyocytes","Spinal_cord","Paraxial_mesoderm","Gut")
celltypes.to.plot <- c("Forebrain_Midbrain_Hindbrain","Neural_crest","Spinal_cord","Erythroid")

diff_filt.dt <- diff.dt[gene%in%genes.to.plot & celltype%in%celltypes.to.plot]

to.plot <- expand.grid(X = unique(diff_filt.dt$celltype), Y = opts$ko.classes) %>% 
  as.data.table %>% setnames(c("celltype","class")) %>%
  merge(diff_filt.dt, by=c("celltype","class"), all.x=T) %>%
  .[complete.cases(.)] %>% setnames("diff","logFC") %>% 
  .[,celltype:=factor(celltype,levels=celltypes.to.plot)] %>% .[,class:=factor(class,levels=opts$ko.classes)]

logFC.limits <- c(-1,4.5)
to.plot %>% .[logFC<logFC.limits[1],logFC:=logFC.limits[1]] %>% .[logFC>logFC.limits[2],logFC:=logFC.limits[2]]

# to.plot.n <-  to.plot %>% 
#   merge(celltype_numbers.dt, by=c("celltype","class")) %>% 
#   .[gene==genes.to.plot[1]] %>%
#   .[,celltype:=factor(celltype,levels=celltypes.to.plot)] %>% .[,class:=factor(class,levels=opts$ko.classes)]

p <- ggplot(to.plot, aes(x=class, y=gene, fill=logFC)) +
  geom_tile(color="black") +
  facet_wrap(~celltype, nrow=1) +
  scale_fill_gradient2(limits = logFC.limits, low = "blue", mid = "white", high = "red", na.value = 'gray70' ) +
  # geom_text(aes(label=Ncells), size=2, data=to.plot.n) +
  geom_text(aes(label=logFC), size=1.75) +
  theme_classic() +
  guides(x = guide_axis(angle = 90)) +
  theme(
    axis.text.x = element_text(color="black", size=rel(0.75)),
    axis.text.y = element_text(color="black", size=rel(0.9)),
    axis.title = element_blank(),
    strip.background = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )


pdf(file.path(io$outdir,"ExE_genes_heatmap_diff.pdf"), width=8, height=3)
print(p)
dev.off()

###############
## HOX genes ##
###############

genes.to.plot <- c("Hoxc9","Hoxc8","Hoxb9","Hoxa9")
# celltypes.to.plot <- c("ExE_mesoderm", "Caudal_Mesoderm", "NMP", "Somitic_mesoderm")
celltypes.to.plot <- c("NMP", "Somitic_mesoderm","Intermediate_mesoderm","ExE_mesoderm")

diff_filt.dt <- diff.dt[gene%in%genes.to.plot & celltype%in%celltypes.to.plot]

to.plot <- expand.grid(X = unique(diff_filt.dt$celltype), Y = opts$ko.classes) %>% 
  as.data.table %>% setnames(c("celltype","class")) %>%
  merge(diff_filt.dt, by=c("celltype","class"), all.x=T) %>%
  .[complete.cases(.)] %>% setnames("diff","logFC") %>% 
  .[,celltype:=factor(celltype,levels=celltypes.to.plot)] %>% .[,class:=factor(class,levels=opts$ko.classes)]

logFC.limits <- c(-2.75,1)
to.plot %>% .[logFC<logFC.limits[1],logFC:=logFC.limits[1]] %>% .[logFC>logFC.limits[2],logFC:=logFC.limits[2]]

# to.plot.n <-  to.plot %>% 
#   merge(celltype_numbers.dt, by=c("celltype","class")) %>% 
#   .[gene==genes.to.plot[1]] %>%
#   .[,celltype:=factor(celltype,levels=celltypes.to.plot)] %>% .[,class:=factor(class,levels=opts$ko.classes)]

p <- ggplot(to.plot, aes(x=class, y=gene, fill=logFC)) +
  geom_tile(color="black") +
  facet_wrap(~celltype, nrow=1) +
  scale_fill_gradient2(limits = logFC.limits, low = "blue", mid = "white", high = "red", na.value = 'gray70' ) +
  # geom_text(aes(label=Ncells), size=2, data=to.plot.n) +
  geom_text(aes(label=logFC), size=1.75) +
  theme_classic() +
  guides(x = guide_axis(angle = 90)) +
  theme(
    axis.text.x = element_text(color="black", size=rel(0.75)),
    axis.text.y = element_text(color="black", size=rel(0.9)),
    axis.title = element_blank(),
    strip.background = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

pdf(file.path(io$outdir,"Hox_genes_heatmap_diff.pdf"), width=8, height=3)
print(p)
dev.off()

######################################
## Heatmap naive pluripotency genes ##
######################################

genes.to.plot <-  c("Zfp42", "Dppa5a", "Dppa3", "Dppa4", "Tfap2c")
# celltypes.to.plot <- c( "Gut", "Cardiomyocytes", "Paraxial_mesoderm", "Spinal_cord")
celltypes.to.plot <- c("Forebrain_Midbrain_Hindbrain", "Spinal_cord","Mesenchyme","Cardiomyocytes")

diff_filt.dt <- diff.dt[gene%in%genes.to.plot & celltype%in%celltypes.to.plot]

to.plot <- expand.grid(X = unique(diff_filt.dt$celltype), Y = opts$ko.classes) %>% 
  as.data.table %>% setnames(c("celltype","class")) %>%
  merge(diff_filt.dt, by=c("celltype","class"), all.x=T) %>%
  .[complete.cases(.)] %>% setnames("diff","logFC") %>% 
  .[,celltype:=factor(celltype,levels=opts$celltypes)] %>% .[,class:=factor(class,levels=opts$ko.classes)]

logFC.limits <- c(-1,2.5)
to.plot %>% .[logFC<logFC.limits[1],logFC:=logFC.limits[1]] %>% .[logFC>logFC.limits[2],logFC:=logFC.limits[2]]

# to.plot.n <-  to.plot %>% 
#   merge(celltype_numbers.dt, by=c("celltype","class")) %>% 
#   .[gene==genes.to.plot[1]] %>%
#   .[,celltype:=factor(celltype,levels=celltypes.to.plot)] %>% .[,class:=factor(class,levels=opts$ko.classes)]

p <- ggplot(to.plot, aes(x=class, y=gene, fill=logFC)) +
  geom_tile(color="black") +
  facet_wrap(~celltype, nrow=1) +
  scale_fill_gradient2(limits = logFC.limits, low = "blue", mid = "white", high = "red", na.value = 'gray70' ) +
  theme_classic() +
  guides(x = guide_axis(angle = 90)) +
  # geom_text(aes(label=N), size=2, data=to.plot.n) +
  geom_text(aes(label=logFC), size=1.75) +
  theme(
    axis.text.x = element_text(color="black", size=rel(0.75)),
    axis.text.y = element_text(color="black", size=rel(0.9)),
    axis.title = element_blank(),
    strip.background = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

pdf(file.path(io$outdir,"Naive_pluripotency_heatmap_diff.pdf"), width=8, height=3)
print(p)
dev.off()

#######################################
## Heatmap primed pluripotency genes ##
#######################################

# Nanog/Sox2 ??

genes.to.plot <-  c("Pou5f1", "Utf1", "Pim2", "Slc7a3", "Fgf5","Gng3","Nanog") # "Dnmt3b"
# celltypes.to.plot <- c( "Gut", "Cardiomyocytes", "Paraxial_mesoderm", "Spinal_cord")
# celltypes.to.plot <- c("Forebrain_Midbrain_Hindbrain", "Spinal_cord","Mesenchyme","Cardiomyocytes")
celltypes.to.plot <- c("Forebrain_Midbrain_Hindbrain","Spinal_cord","NMP","Gut")

diff_filt.dt <- diff.dt[gene%in%genes.to.plot & celltype%in%celltypes.to.plot]

to.plot <- expand.grid(X = unique(diff_filt.dt$celltype), Y = opts$ko.classes) %>% 
  as.data.table %>% setnames(c("celltype","class")) %>%
  merge(diff_filt.dt, by=c("celltype","class"), all.x=T) %>%
  .[complete.cases(.)] %>% setnames("diff","logFC") %>% 
  .[,celltype:=factor(celltype,levels=celltypes.to.plot)] %>% .[,class:=factor(class,levels=opts$ko.classes)]

logFC.limits <- c(-2.5,2.5)
to.plot %>% .[logFC<logFC.limits[1],logFC:=logFC.limits[1]] %>% .[logFC>logFC.limits[2],logFC:=logFC.limits[2]]

# to.plot.n <-  to.plot %>% 
#   merge(celltype_numbers.dt, by=c("celltype","class")) %>% 
#   .[gene==genes.to.plot[1]] %>%
#   .[,celltype:=factor(celltype,levels=celltypes.to.plot)] %>% .[,class:=factor(class,levels=opts$ko.classes)]

p <- ggplot(to.plot, aes(x=class, y=gene, fill=logFC)) +
  geom_tile(color="black") +
  facet_wrap(~celltype, nrow=1) +
  scale_fill_gradient2(limits = logFC.limits, low = "blue", mid = "white", high = "red", na.value = 'gray70' ) +
  theme_classic() +
  guides(x = guide_axis(angle = 90)) +
  # geom_text(aes(label=Ncells), size=2, data=to.plot.n) +
  geom_text(aes(label=logFC), size=1.75) +
  theme(
    axis.text.x = element_text(color="black", size=rel(0.75)),
    axis.text.y = element_text(color="black", size=rel(0.9)),
    axis.title = element_blank(),
    strip.background = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

pdf(file.path(io$outdir,"Primed_pluripotency_heatmap_diff.pdf"), width=8, height=3)
print(p)
dev.off()
