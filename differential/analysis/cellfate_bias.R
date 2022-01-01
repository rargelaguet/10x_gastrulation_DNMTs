here::i_am("differential/analysis/barplots.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

# Load utils
source(here::here("differential/analysis/utils.R"))

##############
## Settings ##
##############

# I/O
io$indir <- file.path(io$basedir,"results_new/differential")
io$outdir <- file.path(io$basedir,"results_new/differential/pdf/cellfate_bias"); dir.create(io$outdir, showWarnings = F)

# Options
opts$min.cells <- 50

###############
## Load data ##
###############

source(here::here("differential/analysis/load_data.R"))

#######################
## Load marker genes ##
#######################

opts$min.marker.score <- 0.85

marker_genes.dt <- fread(io$atlas.marker_genes) %>%
  .[score>=opts$min.marker.score] %>% 
  .[,ExE:=celltype%in%opts$ExE.celltypes]

####################
## Filter results ##
####################

# Filter
diff.dt <- diff.dt[groupA_N>=opts$min.cells & groupB_N>=opts$min.cells]

# Remove sex-specific differences
# dt.filt <- dt[gene!="Xist"]

# Remove Riks
# dt.filt <- dt.filt[!grep("Rik",gene)]

# Remove hemoglobins
# dt.filt <- dt.filt[!grep("^Hb",gene)]

# Subset to lineage markers
diff_markers.dt <- diff.dt[gene%in%unique(marker_genes.dt$gene)]

# Filter out genes that are DE across all celltypes
diff_markers.dt <- diff_markers.dt[,N:=sum(sig),by=c("gene","class")] %>% .[N<=10] %>% .[,N:=NULL]
# diff_markers.dt[class=="E8.5_Dnmt1KO",sum(sig),by=c("gene","class")] %>% View
# diff_markers.dt[class=="E8.5_Dnmt1KO" & gene=="Hoxc8"]

######################################################
## Quantify cell fate disruption using the DE genes ##
######################################################

to.plot <- diff_markers.dt %>%
  merge(
    marker_genes.dt[,c("gene","celltype")] %>% setnames("celltype","celltype_marker"), by = "gene", allow.cartesian=TRUE
  ) %>% .[,sum(sig), by=c("celltype","celltype_marker","class")]

to.plot[class=="E8.5_Dnmt1KO"] %>% View

p <- ggplot(to.plot, aes(x=celltype, y=V1)) +
  geom_bar(aes(fill = celltype_marker), color="black", stat="identity") + 
  facet_wrap(~class, scales="free_y") +
  scale_fill_manual(values=opts$celltype.colors, drop=F) +
  labs(x="", y="Number of DE genes") +
  guides(x = guide_axis(angle = 90)) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.line = element_blank(),
    axis.text.x = element_text(color="black", size=rel(0.65))
  )

pdf(sprintf("%s/DE_barplots_marker_genes_fate_bias.pdf",io$outdir), width=9, height=5)
print(p)
dev.off()

###################################################################
## Quantify cell fate disruption using the DE genes: polar plots ##
###################################################################

for (i in opts$ko.classes) {

  to.plot <- diff_markers.dt[class==i] %>%
    merge(
      marker_genes.dt[,c("gene","celltype")] %>% setnames("celltype","celltype_marker"), by = "gene", allow.cartesian=TRUE
    ) %>% .[,sum(sig), by=c("celltype","celltype_marker","class")]
  
  to.plot[V1>50,V1:=50]
  
  p <- ggplot(to.plot, aes(x=celltype_marker, y=V1)) +
    geom_bar(aes(fill = celltype_marker), color="black", stat = 'identity') + 
    # geom_polygon(color="black", fill=NA, alpha=0.5, linetype="dashed", data=foo) +
    facet_wrap(~celltype) +
    scale_fill_manual(values=opts$celltype.colors, drop=F) +
    coord_polar() + ylim(c(0,50)) +
    # guides(colour = guide_legend(override.aes = list(size=2), ncol=1)) +
    theme_bw() +
    theme(
      legend.position = "none",
      legend.text = element_text(size=rel(0.75)),
      legend.title = element_blank(),
      axis.title=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      axis.line=element_blank(),
      axis.text.x = element_blank()
      # axis.text.x = element_text(angle= -76 - 360 / length(unique(to.plot$celltype)) * seq_along(to.plot$celltype))
    )
  
  pdf(file.path(io$outdir,sprintf("%s_DE_polar_plot_marker_genes_fate_bias.pdf",i)), width=11, height=7)
  print(p)
  dev.off()
  
}

################################################
## Quantify fate disruption: embryonic vs ExE ##
################################################

# Remove genes that are markers of both embryonic and extraembryonic cell types
marker_genes_ExE.dt <- marker_genes.dt %>% 
  .[,V1:=mean(ExE),by="gene"] %>% .[V1%in%c(0,1)] %>% .[,V1:=NULL]

to.plot <- diff_markers.dt %>%
  merge(marker_genes_ExE.dt[,c("gene","ExE")], by="gene", allow.cartesian=TRUE) %>%
  .[,sum(sig), by=c("ExE","celltype","class")]

p <- ggplot(to.plot, aes(x=celltype, y=V1)) +
  geom_bar(aes(fill = ExE), color="black", stat="identity") + 
  facet_wrap(~class, scales="fixed") +
  # scale_fill_manual(values=opts$celltype.colors, drop=F) +
  labs(x="", y="Number of DE genes") +
  guides(x = guide_axis(angle = 90)) +
  theme_bw() +
  theme(
    legend.position = "top",
    axis.line = element_blank(),
    axis.text.x = element_text(color="black", size=rel(0.65))
  )

pdf(sprintf("%s/DE_barplots_marker_genes_fate_bias_ExE.pdf",io$outdir), width=9, height=5.5)
print(p)
dev.off()

