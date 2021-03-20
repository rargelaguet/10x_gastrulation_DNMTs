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
  # "Caudal_epiblast",
  # "PGC",
  # "Anterior_Primitive_Streak",
  # "Notochord",
  # "Def._endoderm",
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
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm"
  # "Parietal_endoderm"
)

# Define groups
opts$groupA <- c(
  # "E8.5_Dnmt3aKO_Dnmt3bWT", 
  "E8.5_Dnmt3aHET_Dnmt3bKO",
  # "E8.5_Dnmt3aHET_Dnmt3bWT", 
  # "E8.5_Dnmt3aKO_Dnmt3bHET", 
  # "E8.5_Dnmt3aKO_Dnmt3bKO"
  "E8.5_Dnmt3aWT_Dnmt3bKO"
)
# "E8.5_Dnmt3aWT_Dnmt3bWT",

opts$groupB <- c("E8.5_Dnmt3aWT_Dnmt3bWT")

#############################################
## Load results from differential analysis ##
#############################################

dt <- opts$groupA %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- sprintf("%s/%s_vs_%s_%s.txt.gz", io$diff.dir,i,opts$groupB,j)
  if (file.exists(file)) fread(file) %>% .[,c(1,2,4,10,11)] %>% .[,c("celltype","groupA","groupB"):=list(j,i,opts$groupB)]
}) %>% rbindlist }) %>% rbindlist

unique(dt$groupA)
unique(dt$groupB)
unique(dt$celltype)

dt[celltype=="Neural_crest" & sig==TRUE] %>% .[,N:=.N,by="gene"] %>% .[N>1] %>% View
####################
## Filter results ##
####################

# Remove some genes
# dt <- dt[gene!="Xist"]

################
## Polar plot ##
################

to.plot <- dt %>% copy %>%
  .[,.(N=sum(sig,na.rm=T)) ,by=c("celltype","groupA","groupB")]

ggplot(to.plot, aes(x=factor(celltype), y=N)) +
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

pdf(sprintf("%s/DE_polar_plots.pdf",io$outdir), width=9, height=5)
print(p)
dev.off()

##############
## Bar plot ##
##############

to.plot <- dt %>% copy %>%
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
    axis.text.x = element_text(color="black", size=rel(0.75))
  )


pdf(sprintf("%s/DE_barplots.pdf",io$outdir), width=9, height=5)
print(p)
dev.off()