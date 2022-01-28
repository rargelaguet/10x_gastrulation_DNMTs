source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$indir <- file.path(io$basedir,"results_all/differential/pseudobulk")
io$outdir <- file.path(io$basedir,"results_all/differential/pseudobulk/pdf/fig"); dir.create(io$outdir, showWarnings = F)

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


##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  setnames("celltype.mapped","celltype") %>%
  .[pass_rnaQC==TRUE & celltype%in%opts$celltypes & class%in%opts$ko.classes] %>%
  .[,celltype:=stringr::str_replace_all(celltype,opts$rename_celltypes)]

celltype_numbers.dt <- sample_metadata[,.N,by=c("class","celltype")]

##############################
## Load precomputed results ##
##############################

# i <- opts$ko.classes[1]; j <- opts$celltypes[10]
diff.dt <- opts$ko.classes %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- sprintf("%s/%s/%s_%s_vs_%s.txt.gz", io$indir,i,j,opts$wt.class,i)
  if (file.exists(file)) {
    fread(file) %>% .[,c("celltype","class"):=list(j,i)]
  }
}) %>% rbindlist }) %>% rbindlist# %>%
  # .[,celltype:=factor(celltype,levels=opts$celltypes)]# %>%
  # .[,class:=factor(class,levels=rev(opts$ko.classes))]

# Merge celltypes 
diff.dt <- diff.dt %>% 
  .[,celltype:=stringr::str_replace_all(celltype,opts$rename_celltypes)] %>%
  .[,.(expr_ko=mean(expr_ko), expr_wt=mean(expr_wt), diff=mean(diff)), by=c("celltype","class","gene")]
  
# save
# fwrite(diff.dt, file.path(io$outdir,"diff_rna_pseudobulk.txt.gz"))

######################
## Heatmap per gene ##
######################

# genes.to.plot <- c("Pim2", "Pou5f1", "Slc7a3", "Utf1", "Dppa5a")
# 
# # i <- "Hoxc8"
# for (i in genes.to.plot) {
#   
#   to.plot <- expand.grid(X = unique(diff.dt$celltype), Y = unique(diff.dt$class)) %>% 
#     as.data.table %>% setnames(c("celltype","class")) %>%
#     merge(diff.dt[gene==i], by=c("celltype","class"), all.x=T)
#   
#   p <- ggplot(to.plot, aes(x=celltype, y=class, fill=diff)) +
#     geom_tile(color="black") +
#     # scale_fill_gradientn(colours = terrain.colors(10), na.value = 'gray70') +
#     # scale_fill_manual(values = wes_palette("GrandBudapest1", n = 3))
#     scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = 'gray70' ) +
#     theme_classic() +
#     guides(x = guide_axis(angle = 90)) +
#     theme(
#       axis.text = element_text(color="black", size=rel(0.7)),
#       axis.title = element_blank(),
#       strip.background = element_blank(),
#       axis.ticks = element_blank(),
#       axis.line = element_blank(),
#       legend.title = element_blank()
#     )
#   
#   
#   pdf(file.path(io$outdir,sprintf("%s.pdf",i)), width=8, height=5)
#   print(p)
#   dev.off()
# }

#######################
## Heatmap ExE genes ##
#######################

genes.to.plot <- c("Rhox9","Trap1a","Xlr3a","Ascl2","Apoe","Fmr1nb","Tex19.1") # "Trh","Rhox5","Rhox6"
# celltypes.to.plot <- c("Caudal_epiblast", "Gut", "Cardiomyocytes", "Rostral_neurectoderm", "Paraxial_mesoderm", "Spinal_cord")
celltypes.to.plot <- c( "Gut", "Cardiomyocytes", "Paraxial_mesoderm", "Spinal_cord")

diff_filt.dt <- diff.dt[gene%in%genes.to.plot & celltype%in%celltypes.to.plot]

to.plot <- expand.grid(X = unique(diff_filt.dt$celltype), Y = opts$ko.classes) %>% 
  as.data.table %>% setnames(c("celltype","class")) %>%
  merge(diff_filt.dt, by=c("celltype","class"), all.x=T) %>%
  .[complete.cases(.)] %>% setnames("diff","logFC") %>% 
  .[,celltype:=factor(celltype,levels=opts$celltypes)] %>% .[,class:=factor(class,levels=opts$ko.classes)]

logFC.limits <- c(-1,5)
to.plot %>% .[logFC<logFC.limits[1],logFC:=logFC.limits[1]] %>% .[logFC>logFC.limits[2],logFC:=logFC.limits[2]]

to.plot.n <-  to.plot %>% merge(celltype_numbers.dt, by=c("celltype","class")) %>% .[gene==genes.to.plot[1]]

p <- ggplot(to.plot, aes(x=class, y=gene, fill=logFC)) +
  geom_tile(color="black") +
  facet_wrap(~celltype, nrow=1) +
  scale_fill_gradient2(limits = logFC.limits, low = "blue", mid = "white", high = "red", na.value = 'gray70' ) +
  geom_text(aes(label=N), size=2, data=to.plot.n) +
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
celltypes.to.plot <- c("ExE_mesoderm", "Caudal_Mesoderm", "NMP", "Somitic_mesoderm")

diff_filt.dt <- diff.dt[gene%in%genes.to.plot & celltype%in%celltypes.to.plot]

to.plot <- expand.grid(X = unique(diff_filt.dt$celltype), Y = opts$ko.classes) %>% 
  as.data.table %>% setnames(c("celltype","class")) %>%
  merge(diff_filt.dt, by=c("celltype","class"), all.x=T) %>%
  .[complete.cases(.)] %>% setnames("diff","logFC") %>% 
  .[,celltype:=factor(celltype,levels=opts$celltypes)] %>% .[,class:=factor(class,levels=opts$ko.classes)]

logFC.limits <- c(-2.75,1)
to.plot %>% .[logFC<logFC.limits[1],logFC:=logFC.limits[1]] %>% .[logFC>logFC.limits[2],logFC:=logFC.limits[2]]

to.plot.n <-  to.plot %>% merge(celltype_numbers.dt, by=c("celltype","class")) %>% .[gene==genes.to.plot[1]]

p <- ggplot(to.plot, aes(x=class, y=gene, fill=logFC)) +
  geom_tile(color="black") +
  facet_wrap(~celltype, nrow=1) +
  scale_fill_gradient2(limits = logFC.limits, low = "blue", mid = "white", high = "red", na.value = 'gray70' ) +
  geom_text(aes(label=N), size=2, data=to.plot.n) +
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

genes.to.plot <-  c("Zfp42", "Dppa5a", "Dppa3", "Dppa4", "Tfap2c", "Pecam1")
celltypes.to.plot <- c( "Gut", "Cardiomyocytes", "Paraxial_mesoderm", "Spinal_cord")

diff_filt.dt <- diff.dt[gene%in%genes.to.plot & celltype%in%celltypes.to.plot]

to.plot <- expand.grid(X = unique(diff_filt.dt$celltype), Y = opts$ko.classes) %>% 
  as.data.table %>% setnames(c("celltype","class")) %>%
  merge(diff_filt.dt, by=c("celltype","class"), all.x=T) %>%
  .[complete.cases(.)] %>% setnames("diff","logFC") %>% 
  .[,celltype:=factor(celltype,levels=opts$celltypes)] %>% .[,class:=factor(class,levels=opts$ko.classes)]

logFC.limits <- c(-1,2.5)
to.plot %>% .[logFC<logFC.limits[1],logFC:=logFC.limits[1]] %>% .[logFC>logFC.limits[2],logFC:=logFC.limits[2]]

to.plot.n <-  to.plot %>% merge(celltype_numbers.dt, by=c("celltype","class")) %>% .[gene==genes.to.plot[1]]

p <- ggplot(to.plot, aes(x=class, y=gene, fill=logFC)) +
  geom_tile(color="black") +
  facet_wrap(~celltype, nrow=1) +
  scale_fill_gradient2(limits = logFC.limits, low = "blue", mid = "white", high = "red", na.value = 'gray70' ) +
  theme_classic() +
  guides(x = guide_axis(angle = 90)) +
  geom_text(aes(label=N), size=2, data=to.plot.n) +
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

genes.to.plot <-  c("Pou5f1", "Utf1", "Pim2", "Slc7a3", "Fgf5","Gng3") # "Dnmt3b"
celltypes.to.plot <- c( "Gut", "Cardiomyocytes", "Paraxial_mesoderm", "Spinal_cord")

diff_filt.dt <- diff.dt[gene%in%genes.to.plot & celltype%in%celltypes.to.plot]

to.plot <- expand.grid(X = unique(diff_filt.dt$celltype), Y = opts$ko.classes) %>% 
  as.data.table %>% setnames(c("celltype","class")) %>%
  merge(diff_filt.dt, by=c("celltype","class"), all.x=T) %>%
  .[complete.cases(.)] %>% setnames("diff","logFC") %>% 
  .[,celltype:=factor(celltype,levels=opts$celltypes)] %>% .[,class:=factor(class,levels=opts$ko.classes)]

logFC.limits <- c(-1,2.25)
to.plot %>% .[logFC<logFC.limits[1],logFC:=logFC.limits[1]] %>% .[logFC>logFC.limits[2],logFC:=logFC.limits[2]]

to.plot.n <-  to.plot %>% merge(celltype_numbers.dt, by=c("celltype","class")) %>% .[gene==genes.to.plot[1]]

p <- ggplot(to.plot, aes(x=class, y=gene, fill=logFC)) +
  geom_tile(color="black") +
  facet_wrap(~celltype, nrow=1) +
  scale_fill_gradient2(limits = logFC.limits, low = "blue", mid = "white", high = "red", na.value = 'gray70' ) +
  theme_classic() +
  guides(x = guide_axis(angle = 90)) +
  geom_text(aes(label=N), size=2, data=to.plot.n) +
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
