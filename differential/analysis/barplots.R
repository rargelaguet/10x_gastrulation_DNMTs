##############
## Settings ##
##############

# Load settings
source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
source("/Users/ricard/10x_gastrulation_DNMTs/differential/analysis/utils.R")

# I/O
io$diff.dir <- paste0(io$basedir,"/results/differential")
io$outdir <- paste0(io$basedir,"/results/differential/pdf")

opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  "PGC",
  "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  "Mixed_mesoderm",
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
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm"
  # "ExE_endoderm",
  # "ExE_ectoderm",
  # "Parietal_endoderm"
)

opts$ExE.celltypes <- c(
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm",
  "Parietal_endoderm"
)

# Define groups
opts$groupA <- c(
  "E8.5_Dnmt3aKO_Dnmt3bWT",
  "E8.5_Dnmt3aHET_Dnmt3bKO",
  # "E8.5_Dnmt3aHET_Dnmt3bWT",
  "E8.5_Dnmt3aKO_Dnmt3bHET",
  "E8.5_Dnmt3aKO_Dnmt3bKO",
  "E8.5_Dnmt3aWT_Dnmt3bKO",
  "E8.5_Dnmt1KO"
)

opts$groupB <- c("E8.5_WT")

opts$min.cells <- 50

#############################################
## Load results from differential analysis ##
#############################################

dt <- opts$groupA %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- sprintf("%s/%s_vs_%s_%s.txt.gz", io$diff.dir,i,opts$groupB,j)
  if (file.exists(file)) fread(file) %>% .[,c(1,2,4,6,7,10,11)] %>% .[,c("celltype","groupA","groupB"):=list(j,i,opts$groupB)]
}) %>% rbindlist }) %>% rbindlist

dt[,celltype:=factor(celltype,levels=opts$celltypes)]
dt[,groupA:=factor(groupA,levels=opts$groupA)]

unique(dt$groupA)
unique(dt$groupB)
unique(dt$celltype)

# dt[celltype=="Neural_crest" & sig==TRUE] %>% .[,N:=.N,by="gene"] %>% .[N>1] %>% View
# dt.filt[groupA=="E8.5_Dnmt3aKO_Dnmt3bWT" & celltype=="Caudal_epiblast" & sig==TRUE] %>% View

####################
## Filter results ##
####################

# Filter
dt.filt <- dt[groupA_N>=opts$min.cells & groupB_N>=opts$min.cells]

# Remove sex-specific differences
# dt.filt <- dt[gene!="Xist"]

# Remove Riks
dt.filt <- dt.filt[!grep("Rik",gene)]

# Remove hemoglobins
dt.filt <- dt.filt[!grep("^Hb",gene)]

# Subset to lineage markers
marker_genes.dt <- fread(io$atlas.marker_genes) %>% .[,ExE:=celltype%in%opts$ExE.celltypes]
dt.filt <- dt.filt[gene%in%unique(marker_genes.dt$gene)]
# marker_genes.dt[gene=="H19"]

################
## Polar plot ##
################

to.plot <- dt.filt %>% copy %>%
  .[,.(N=sum(sig,na.rm=T)) ,by=c("celltype","groupA","groupB")]

p <- ggplot(to.plot, aes(x=factor(celltype), y=N)) +
  geom_bar(aes(fill = celltype), color="black", stat = 'identity') + 
  # geom_polygon(color="black", fill=NA, alpha=0.5, linetype="dashed", data=foo) +
  facet_wrap(~groupA) +
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

pdf(sprintf("%s/DE_polar_plots_marker_genes.pdf",io$outdir), width=9, height=5)
print(p)
dev.off()

##############
## Bar plot ##
##############

to.plot <- dt.filt %>% copy %>%
  .[,.(N=sum(sig,na.rm=T)) ,by=c("celltype","groupA","groupB")]

p <- ggplot(to.plot, aes(x=factor(celltype), y=N)) +
  geom_bar(aes(fill = celltype), color="black", stat = 'identity') + 
  facet_wrap(~groupA) +
  scale_fill_manual(values=opts$celltype.colors, drop=F) +
  labs(x="", y="Number of DE genes") +
  guides(x = guide_axis(angle = 90)) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.line = element_blank(),
    axis.text.x = element_text(color="black", size=rel(0.65))
  )


pdf(sprintf("%s/DE_barplots_marker_genes.pdf",io$outdir), width=9, height=5)
print(p)
dev.off()

##################################
## Quantify fate disruption (?) ##
##################################

to.plot <- dt.filt[sig==T] %>%
  merge(
    marker_genes.dt[,c("gene","celltype","ExE")] %>% setnames("celltype","celltype_marker"),
    by = c("gene"), allow.cartesian=TRUE
  ) %>%
  .[,.N, by=c("celltype","celltype_marker","groupA","groupB")]

p <- ggplot(to.plot, aes(x=factor(celltype), y=N)) +
  geom_bar(aes(fill = celltype_marker), color="black", stat="identity") + 
  facet_wrap(~groupA, scales="free_y") +
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

################################################
## Quantify fate disruption: embryonic vs ExE ##
################################################

to.plot <- dt.filt[sig==T] %>%
  merge(
    marker_genes.dt[,c("gene","ExE")],
    by = c("gene"), allow.cartesian=TRUE
  ) %>%
  .[,.N, by=c("ExE","celltype","groupA","groupB")]

p <- ggplot(to.plot, aes(x=factor(celltype), y=N)) +
  geom_bar(aes(fill = ExE), color="black", stat="identity") + 
  facet_wrap(~groupA, scales="free_y") +
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

##########################################
## Barplot of fraction of genes up/down ##
##########################################

stop("DOUBLE CHECK THE DIRECTION")
to.plot <- dt.filt[sig==T] %>%
  .[,direction:=c("Down","Up")[as.numeric(logFC>0)+1]] %>%
  .[,.N, by=c("celltype","groupA","groupB","direction")]

p <- ggplot(to.plot, aes(x=factor(celltype), y=N)) +
  geom_bar(aes(fill = direction), color="black", stat="identity") + 
  facet_wrap(~groupA, scales="free_y") +
  labs(x="", y="Number of DE genes") +
  guides(x = guide_axis(angle = 90)) +
  theme_bw() +
  theme(
    legend.position = "top",
    axis.line = element_blank(),
    axis.text.x = element_text(color="black", size=rel(0.65))
  )

pdf(sprintf("%s/DE_barplots_marker_genes_direction.pdf",io$outdir), width=9, height=5.5)
print(p)
dev.off()

################################################
## Identify genes with unidirectional effects ##
################################################

stop("DOUBLE CHECK THE DIRECTION")

foo <- dt.filt[sig==T] %>%
  .[,direction:=c("Down","Up")[as.numeric(logFC>0)+1]] %>%
  .[,.N, by=c("gene","groupA","groupB","direction")] %>%
  dcast(gene+groupA+groupB~direction,value.var="N", fill=0)

# Plot number of Up/Down/Mixed genes per class
foo[,class:="Down"]
foo[Up>0 & Down==0, class:="Up"]
foo[Up>0 & Down>0, class:="Mixed"]

to.plot <- foo[,.N,by=c("class","groupA")]

p <- ggplot(foo[,.N,by=c("class","groupA")], aes(x=class, y=N)) +
  geom_bar(fill="gray70", color="black", stat="identity") + 
  facet_wrap(~groupA, scales="free_y") +
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
to.plot <- foo %>% 
  split(.$groupA) %>% map(~ setorder(., -Up) %>% head(n=30)) %>% rbindlist  
# .[,gene:=factor(gene,levels=gene)]

p <- ggplot(to.plot, aes_string(x="gene", y="Up")) +
  facet_wrap(~groupA, scales="free_y") +
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
to.plot <- foo %>% 
  split(.$groupA) %>% map(~ setorder(., -Down) %>% head(n=30)) %>% rbindlist  
# .[,gene:=factor(gene,levels=gene)]

p <- ggplot(to.plot, aes_string(x="gene", y="Down")) +
  facet_wrap(~groupA, scales="free_y") +
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

#############
## Explore ##
#############

foo <- dt.filt[sig==T] %>%
  merge(
    marker_genes.dt[,c("gene","celltype")] %>% setnames("celltype","celltype_marker"),
    by = c("gene"), allow.cartesian=TRUE
  )

bar <- foo[celltype_marker=="Visceral_endoderm"] %>% 
  .[,.N,by=c("gene","groupA","groupB")]


dt.filt[gene=="Rhox5" & sig==TRUE] %>% View
