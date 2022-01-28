here::i_am("differential/analysis/barplots.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

##############
## Settings ##
##############

# I/O
io$indir <- file.path(io$basedir,"results_all/differential")
io$outdir <- file.path(io$basedir,"results_all/differential/pdf/fig"); dir.create(io$outdir, showWarnings = F)

# Options
opts$min.cells <- 50

###############
## Load data ##
###############

source(here::here("differential/analysis/load_data.R"))

diff.dt <- diff.dt[!is.na(logFC)]

####################
## Filter results ##
####################

# Filter
diff.dt <- diff.dt[groupA_N>=opts$min.cells & groupB_N>=opts$min.cells]

# Remove sex-specific differences
diff.dt <- diff.dt[gene!="Xist"]

# Filter out genes manually
# dt.filt <- dt.filt[!grep("[^Rik|^Gm|^Rpl]",gene)]
diff.dt <- diff.dt[!grepl("^Hb",gene)]
diff.dt <- diff.dt[gene!="Cdkn1c"]

# Subset to lineage markers
marker_genes.dt <- fread(io$atlas.marker_genes) %>% .[,ExE:=celltype%in%opts$ExE.celltypes]
diff_markers.dt <- diff.dt[gene%in%unique(marker_genes.dt$gene)]

################
## Polar plot ##
################

to.plot <- diff_markers.dt %>% copy %>%
  .[,.(N=sum(sig,na.rm=T)) ,by=c("celltype","class")]

to.plot[N>=80,N:=80] # For viz purposes

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


to.plot <- diff_markers.dt[sig==T] %>%
  .[,direction:=c("Downregulated in KO","Upregulated in KO")[as.numeric(logFC>0)+1]] %>%
  .[,direction:=factor(direction,levels=c("Upregulated in KO","Downregulated in KO"))] %>%
  .[,.N, by=c("celltype","class","direction")]

to.plot[N>=70,N:=70] # For viz purposes

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

marker_genes_ExE.dt <- marker_genes.dt %>% 
  .[,V1:=mean(ExE),by="gene"] %>% .[V1%in%c(0,1)] %>% .[,V1:=NULL]

to.plot <- diff_markers.dt %>%
  .[class=="Dnmt1_KO"] %>%
  merge(marker_genes_ExE.dt[,c("gene","ExE")], by="gene", allow.cartesian=TRUE) %>%
  .[,.(N=sum(sig)), by=c("ExE","celltype","class","sign")] %>%
  .[,ExE:=factor(ExE, levels=c(FALSE,TRUE))]

p <- ggplot(to.plot, aes(x=celltype, y=N)) +
  geom_bar(aes(fill = ExE), color="black", stat="identity") + 
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

pdf(file.path(io$outdir,"DE_barplots_ExE_up_down.pdf"), width=6, height=5)
print(p)
dev.off()

##################################
## Embryonic vs ExE polar plots ##
##################################

to.plot <- diff_markers.dt %>%
  .[class=="Dnmt1_KO" & gene%in%marker_genes_ExE.dt[ExE==TRUE,gene]] %>%
  .[,.(N=sum(sig)), by=c("celltype","sign")]

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

epiblast.markers <- unique(marker_genes.dt[celltype=="Epiblast",gene])

to.plot <- diff_markers.dt %>%
  .[class=="Dnmt1_KO" & gene%in%epiblast.markers] %>%
  .[,.(N=sum(sig)), by=c("celltype","sign")]

p <- ggplot(to.plot, aes(x=celltype, y=N)) +
  geom_bar(aes(fill = celltype), color="black", stat = 'identity') + 
  scale_fill_manual(values=opts$celltype.colors, drop=F) +
  coord_polar() +
  facet_wrap(~sign, scales="fixed", nrow=1) +
  geom_text(aes(label=N), size=5, data=to.plot[N>=5]) +
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
# to.plot <- to.plot[celltype_marker!="Notochord"]

# tmp <- diff_markers.dt[class=="Dnmt1_KO" & sig==T] %>% merge(marker_genes.dt[,c("gene","celltype")] %>% setnames("celltype","celltype_marker"), by = "gene", allow.cartesian=TRUE)

to.plot <- to.plot[V1>=5]
to.plot[,celltype:=factor(celltype, levels=opts$celltypes[opts$celltypes%in%unique(to.plot$celltype)])]

p <- ggplot(to.plot, aes(x=celltype, y=V1)) +
  geom_bar(aes(fill = celltype_marker), color="black", stat="identity") + 
  facet_wrap(~sign, scales="free_y", nrow=2) +
  # facet_wrap(~class, scales="fixed") +
  scale_fill_manual(values=opts$celltype.colors[names(opts$celltype.colors)%in%unique(to.plot$celltype_marker)], drop=F) +
  labs(x="", y="Number of DE genes") +
  guides(x = guide_axis(angle = 90)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    legend.position = "none",
    axis.line = element_blank(),
    axis.text.x = element_text(color="black", size=rel(0.75))
  )

pdf(sprintf("%s/DE_barplots_marker_genes_fate_bias.pdf",io$outdir), width=9, height=5)
print(p)
dev.off()
