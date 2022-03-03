here::i_am("differential/analysis/barplots.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

##############
## Settings ##
##############

# I/O
io$indir <- file.path(io$basedir,"results/differential")
io$outdir <- file.path(io$basedir,"results/differential/pdf/fig"); dir.create(io$outdir, showWarnings = F, recursive = T)

# Options
opts$min.cells <- 50

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

###############
## Load data ##
###############

source(here::here("differential/analysis/load_data.R"))

diff.dt <- diff.dt[!is.na(logFC)]

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  setnames("celltype.mapped","celltype") %>%
  .[,celltype:=stringr::str_replace_all(celltype,opts$rename_celltypes)] %>%
  .[pass_rnaQC==TRUE & celltype%in%opts$celltypes & class%in%c(opts$wt.class,opts$ko.classes)]

celltype_numbers.dt <- sample_metadata[,.(Ncells=.N, Nembryos=length(unique(sample))), by=c("class","celltype")]

stopifnot(unique(diff.dt$celltype)%in%unique(celltype_numbers.dt$celltype))

####################
## Filter results ##
####################

# Filter out all cell types for which we don't have enough observations in the WT
celltypes.to.use <- celltype_numbers.dt[class=="WT" & Nembryos>=3 & Ncells>=150,celltype]
diff.dt <- diff.dt[celltype%in%celltypes.to.use]

# Filter out all cell types for which we don't have enough observations in the KO
# diff.dt <- diff.dt[groupA_N>=opts$min.cells & groupB_N>=opts$min.cells]
diff.dt <- diff.dt %>% merge(celltype_numbers.dt[Nembryos>=3 & Ncells>=50,c("celltype","class")], by=c("class","celltype"))

# Remove sex-specific differences
diff.dt <- diff.dt[gene!="Xist"]

# Filter out genes manually
diff.dt <- diff.dt[!grepl("^Rik",gene)]
diff.dt <- diff.dt[!grepl("Rik$",gene)]
diff.dt <- diff.dt[!grepl("^Hb",gene)]
diff.dt <- diff.dt[gene!="Cdkn1c"]

# Subset to lineage markers
marker_genes.dt <- fread(io$atlas.marker_genes)

# Divide genes into embryonic and ExE
marker_genes_ExE.dt <- marker_genes.dt %>% 
  .[,ExE:=celltype%in%opts$ExE.celltypes] %>%
  .[,V1:=mean(ExE),by="gene"] %>% .[V1<0.30|V1>0.70] %>% .[,ExE:=ifelse(V1>0.5,TRUE,FALSE)]

diff_markers.dt <- diff.dt[gene%in%unique(marker_genes.dt$gene)]

################
## Polar plot ##
################

to.plot <- diff_markers.dt %>% copy %>%
  .[,.(N=sum(sig,na.rm=T)) ,by=c("celltype","class")]

# to.plot[N>=80,N:=80] # For viz purposes

p <- ggplot(to.plot, aes(x=celltype, y=N)) +
  geom_bar(aes(fill = celltype), color="black", stat = 'identity') + 
  facet_wrap(~class) +
  scale_fill_manual(values=opts$celltype.colors, drop=F) +
  coord_polar() +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    legend.position = "none",
    legend.text = element_text(size=rel(0.75)),
    legend.title = element_blank(),
    axis.title=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    axis.line=element_blank(),
    axis.text.x = element_blank()
  )

pdf(file.path(io$outdir,"DE_polar_plots_marker_genes.pdf"), width=11, height=8)
print(p)
dev.off()

#########################################
## Polarplot fraction of genes up/down ##
#########################################

to.plot <- diff_markers.dt %>%
  .[sig==T] %>%
  .[,direction:=c("Downregulated in KO","Upregulated in KO")[as.numeric(logFC>0)+1]] %>%
  .[,direction:=factor(direction,levels=c("Upregulated in KO","Downregulated in KO"))] %>%
  .[,.N, by=c("celltype","class","direction")]

# to.plot[N>=70,N:=70] # For viz purposes

p <- ggplot(to.plot, aes(x=celltype, y=N)) +
  geom_bar(aes(fill = direction), color="black", stat = 'identity') + 
  facet_wrap(~class) +
  coord_polar() +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    legend.position = "top",
    legend.text = element_text(size=rel(0.75)),
    legend.title = element_blank(),
    axis.title=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    axis.line=element_blank(),
    axis.text.x = element_blank()
  )

pdf(file.path(io$outdir,"DE_polar_plots_marker_genes_up_down.pdf"), width=11, height=8)
print(p)
dev.off()

################################
## Embryonic vs ExE bar plots ##
################################

to.plot <- diff_markers.dt %>%
  .[class=="Dnmt1_KO"] %>%
  merge(unique(marker_genes_ExE.dt[,c("gene","ExE")]), by="gene", allow.cartesian=TRUE) %>%
  .[,.(N=sum(sig)), by=c("ExE","celltype","class","sign")] %>%
  .[,ExE:=factor(ExE, levels=c(FALSE,TRUE))]

to.plot <- to.plot[!celltype%in%c(opts$ExE.celltypes)]

p <- ggplot(to.plot, aes(x=celltype, y=N)) +
  geom_bar(aes(fill = ExE), color="black", stat="identity") + 
  facet_wrap(~sign, scales="fixed", nrow=2) +
  # scale_fill_manual(values=opts$celltype.colors, drop=F) +
  labs(x="", y="Number of DE genes") +
  guides(x = guide_axis(angle = 90)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    legend.position = "none",
    axis.line = element_blank(),
    axis.text.x = element_text(color="black", size=rel(0.65)),
    axis.text.y = element_text(color="black", size=rel(0.85))
  )

pdf(file.path(io$outdir,"DE_barplots_ExE_up_down.pdf"), width=6.5, height=4.5)
print(p)
dev.off()

##################################
## Embryonic vs ExE polar plots ##
##################################

to.plot <- diff_markers.dt %>%
  .[class=="Dnmt1_KO" & gene%in%marker_genes_ExE.dt[ExE==TRUE,gene]] %>%
  .[,.(N=sum(sig)), by=c("celltype","sign")]

# to.plot[N>=22,N:=22] # For viz purposes

p <- ggplot(to.plot, aes(x=celltype, y=N)) +
  geom_bar(aes(fill = celltype), color="black", stat = 'identity') +
  scale_fill_manual(values=opts$celltype.colors, drop=F) +
  geom_text(aes(label=N), size=5, data=to.plot[N>=12]) +
  coord_polar() +
  facet_wrap(~sign, scales="fixed", nrow=1) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    legend.position = "none",
    legend.text = element_text(size=rel(0.75)),
    legend.title = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank()
  )

pdf(file.path(io$outdir,"DE_polar_plots_ExE_up_down.pdf"), width=6, height=3)
print(p)
dev.off()


#########################
## Hox genes bar plots ##
#########################

marker_genes_hox.dt <- marker_genes.dt %>% 
  .[,hox:=grepl("^Hox",gene)] %>%
  .[,V1:=mean(hox),by="gene"] %>% .[V1%in%c(0,1)] %>% .[,V1:=NULL]

to.plot <- diff_markers.dt %>%
  .[class=="Dnmt1_KO"] %>%
  merge(marker_genes_hox.dt[,c("gene","hox")], by="gene", allow.cartesian=TRUE) %>%
  .[,.(N=sum(sig)), by=c("hox","celltype","class","sign")] %>%
  .[,hox:=factor(hox, levels=c(FALSE,TRUE))]

p <- ggplot(to.plot, aes(x=celltype, y=N)) +
  geom_bar(aes(fill = hox), color="black", stat="identity") + 
  facet_wrap(~sign, scales="free_y", nrow=2) +
  # scale_fill_manual(values=opts$celltype.colors, drop=F) +
  labs(x="", y="Number of DE genes") +
  guides(x = guide_axis(angle = 90)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    legend.position = "top",
    axis.line = element_blank(),
    axis.text.x = element_text(color="black", size=rel(0.65)),
    axis.text.y = element_text(color="black", size=rel(0.85))
  )

pdf(file.path(io$outdir,"DE_barplots_hox_up_down.pdf"), width=6, height=5)
print(p)
dev.off()

###########################
## Hox genes polar plots ##
###########################

to.plot <- diff_markers.dt %>%
  .[grepl("^Hox",gene) & class=="Dnmt1_KO"] %>%
  .[,.(N=sum(sig)), by=c("celltype","sign")] 

p <- ggplot(to.plot, aes(x=celltype, y=N)) +
  geom_bar(aes(fill = celltype), color="black", stat = 'identity') + 
  scale_fill_manual(values=opts$celltype.colors, drop=F) +
  geom_text(aes(label=N), size=5, data=to.plot[N>=5]) +
  coord_polar() +
  facet_wrap(~sign, scales="fixed", nrow=1) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    legend.position = "none",
    legend.text = element_text(size=rel(0.75)),
    legend.title = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank()
  )

pdf(file.path(io$outdir,"DE_polar_plots_hox_up_down.pdf"), width=6, height=3)
print(p)
dev.off()

#######################################
## Primed pluripotency markers plots ##
#######################################

epiblast.markers <- c(unique(marker_genes.dt[celltype=="Epiblast",gene]),"Nanog")

to.plot <- diff_markers.dt %>%
  .[class=="Dnmt1_KO" & gene%in%epiblast.markers] %>%
  .[,.(N=sum(sig)), by=c("celltype","sign")]

# to.plot[N>=7,N:=7] # For viz purposes

p <- ggplot(to.plot, aes(x=celltype, y=N)) +
  geom_bar(aes(fill = celltype), color="black", stat = 'identity') + 
  scale_fill_manual(values=opts$celltype.colors, drop=F) +
  coord_polar() +
  facet_wrap(~sign, scales="fixed", nrow=1) +
  geom_text(aes(label=N), size=5, data=to.plot[N>=3]) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    legend.position = "none",
    legend.text = element_text(size=rel(0.75)),
    legend.title = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank()
  )

pdf(file.path(io$outdir,"DE_polar_plots_Pluripotency_up_down.pdf"), width=6, height=3)
print(p)
dev.off()

#############################
## Cell fate bias barplots ##
#############################

to.plot <- diff_markers.dt %>%
  .[class=="Dnmt1_KO"] %>%
  merge(
    marker_genes.dt[,c("gene","celltype")] %>% setnames("celltype","celltype_marker"), by = "gene", allow.cartesian=TRUE
  ) %>% .[,sum(sig), by=c("celltype","celltype_marker","class","sign")]

# to.plot <- to.plot[class%in%c("Dnmt1_KO","Dnmt3ab_KO")]
to.plot <- to.plot[celltype_marker!="Notochord"]
to.plot <- to.plot[!celltype%in%c(opts$ExE.celltypes,"Surface_ectoderm")]

# tmp <- diff_markers.dt[class=="Dnmt1_KO" & sig==T] %>% merge(marker_genes.dt[,c("gene","celltype")] %>% setnames("celltype","celltype_marker"), by = "gene", allow.cartesian=TRUE)

# to.plot <- to.plot[V1>=6]
to.plot[,celltype:=factor(celltype, levels=opts$celltypes[opts$celltypes%in%unique(to.plot$celltype)])]

p <- ggplot(to.plot, aes(x=celltype, y=V1)) +
  geom_bar(aes(fill = celltype_marker), color="black", stat="identity") + 
  facet_wrap(~sign, scales="fixed", nrow=2) +
  # facet_wrap(~class, scales="fixed") +
  scale_fill_manual(values=opts$celltype.colors[names(opts$celltype.colors)%in%unique(to.plot$celltype_marker)], drop=F) +
  labs(x="", y="Number of DE genes") +
  guides(x = guide_axis(angle = 90)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    legend.position = "none",
    axis.line = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color="black", size=rel(1)),
    axis.text.x = element_text(color="black", size=rel(0.75))
  )

pdf(sprintf("%s/DE_barplots_marker_genes_fate_bias.pdf",io$outdir), width=8, height=5)
print(p)
dev.off()
