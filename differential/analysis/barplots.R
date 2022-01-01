here::i_am("differential/analysis/barplots.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

# Load utils
# source(here::here("differential/analysis/utils.R"))

##############
## Settings ##
##############

# I/O
io$indir <- file.path(io$basedir,"results_new/differential")
io$outdir <- file.path(io$basedir,"results_new/differential/pdf/volcano_plots"); dir.create(io$outdir, showWarnings = F)

# Options
opts$min.cells <- 50

###############
## Load data ##
###############

source(here::here("differential/analysis/load_data.R"))

####################
## Filter results ##
####################

# Filter
diff.dt <- diff.dt[groupA_N>=opts$min.cells & groupB_N>=opts$min.cells]

# Remove sex-specific differences
# dt.filt <- dt[gene!="Xist"]

# Filter out genes manually
# dt.filt <- dt.filt[!grep("[^Rik|^Gm|^Rpl]",gene)]
# dt.filt <- dt.filt[!grep("^Hb",gene)]

# Subset to lineage markers
marker_genes.dt <- fread(io$atlas.marker_genes)# %>% .[,ExE:=celltype%in%opts$ExE.celltypes]
diff_markers.dt <- diff.dt[gene%in%unique(marker_genes.dt$gene)]

################
## Polar plot ##
################

to.plot <- diff_markers.dt %>% copy %>%
  .[,.(N=sum(sig,na.rm=T)) ,by=c("celltype","class")]

p <- ggplot(to.plot, aes(x=celltype, y=N)) +
  geom_bar(aes(fill = celltype), color="black", stat = 'identity') + 
  # geom_polygon(color="black", fill=NA, alpha=0.5, linetype="dashed", data=foo) +
  facet_wrap(~class) +
  scale_fill_manual(values=opts$celltype.colors, drop=F) +
  coord_polar() +
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

pdf(file.path(io$outdir,"DE_polar_plots_marker_genes.pdf"), width=11, height=8)
print(p)
dev.off()

##############
## Bar plot ##
##############

to.plot <- diff_markers.dt[sig==T] %>%
  .[,.N, by=c("celltype","class")]

p <- ggplot(to.plot, aes(x=celltype, y=N)) +
  geom_bar(aes(fill = celltype), color="black", stat = 'identity') + 
  facet_wrap(~class) +
  scale_fill_manual(values=opts$celltype.colors, drop=F) +
  labs(x="", y="Number of DE genes") +
  guides(x = guide_axis(angle = 90)) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.line = element_blank(),
    axis.text.x = element_text(color="black", size=rel(0.75))
  )


pdf(file.path(io$outdir,"DE_barplots_marker_genes.pdf"), width=11, height=8)
print(p)
dev.off()

##########################################
## Barplot of fraction of genes up/down ##
##########################################

to.plot <- diff_markers.dt[sig==T] %>%
  .[,direction:=c("Downregulated in KO","Upregulated in KO")[as.numeric(logFC>0)+1]] %>%
  .[,.N, by=c("celltype","class","direction")]

p <- ggplot(to.plot, aes(x=factor(celltype), y=N)) +
  geom_bar(aes(fill = direction), color="black", stat="identity") + 
  facet_wrap(~class, scales="free_y") +
  labs(x="", y="Number of DE genes") +
  guides(x = guide_axis(angle = 90)) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_text(color="black", size=rel(0.75))
  )

pdf(file.path(io$outdir,"DE_barplots_marker_genes_direction.pdf"), width=11, height=8)
print(p)
dev.off()

################################################
## Identify genes with unidirectional effects ##
################################################

foo <- diff_markers.dt[sig==T] %>%
  .[,direction:=c("Down","Up")[as.numeric(logFC>0)+1]] %>%
  .[,.N, by=c("gene","class","direction")] %>%
  dcast(gene+class~direction,value.var="N", fill=0)

foo[,diff_type:="Down"]
foo[Up>0 & Down==0, diff_type:="Up"]
foo[Up>0 & Down>0, diff_type:="Mixed"]

to.plot <- foo[,.N,by=c("class","diff_type")]

# Plot (...)
p <- ggplot(to.plot, aes(x=class, y=N, fill=diff_type)) +
  geom_bar(color="black", stat="identity", position="dodge") + 
  # facet_wrap(~class, scales="fixed") +
  labs(x="", y="Number of DE genes") +
  # guides(x = guide_axis(angle = 90)) +
  theme_bw() +
  theme(
    axis.text.y = element_text(color="black", size=rel(1)),
    axis.text.x = element_text(color="black", size=rel(1))
  )

pdf(sprintf("%s/DE_barplots_marker_genes_direction_overall.pdf",io$outdir), width=9, height=5.5)
print(p)
dev.off()


# Plot top DE genes with consistent upregulation
to.plot <- foo %>% split(.$class) %>% map(~ setorder(., -Up) %>% head(n=30)) %>% rbindlist  

p <- ggplot(to.plot, aes_string(x="gene", y="Up")) +
  facet_wrap(~class, scales="free_y") +
  geom_point(size=2) +
  geom_segment(aes_string(xend="gene"), size=0.75, yend=0) +
  coord_flip() +
  labs(y="Number of DE events (Upregulated genes)", x="") +
  theme_bw() +
  theme(
    axis.text = element_text(color="black", size=rel(0.8)),
    axis.ticks = element_line(size=rel(0.6))
  )
  
pdf(sprintf("%s/upregulated_genes.pdf",io$outdir), width=8, height=8)
print(p)
dev.off()

# Plot top DE genes with consistent downregulation
to.plot <- foo %>% split(.$class) %>% map(~ setorder(., -Down) %>% head(n=30)) %>% rbindlist  

p <- ggplot(to.plot, aes_string(x="gene", y="Down")) +
  facet_wrap(~class, scales="free_y") +
  geom_point(size=2) +
  geom_segment(aes_string(xend="gene"), size=0.75, yend=0) +
  coord_flip() +
  labs(y="Number of DE events (Downregulated genes)", x="") +
  theme_bw() +
  theme(
    axis.text = element_text(color="black", size=rel(0.8)),
    axis.ticks = element_line(size=rel(0.6))
  )

pdf(sprintf("%s/downregulated_genes.pdf",io$outdir), width=8, height=8)
print(p)
dev.off()
