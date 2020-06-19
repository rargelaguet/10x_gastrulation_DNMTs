##############
## Settings ##
##############

source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
source("/Users/ricard/10x_gastrulation_DNMTs/differential/analysis/utils.R")
io$diff.dir <- paste0(io$basedir,"/results/differential")
io$outdir <- paste0(io$basedir,"/results/differential/pdf")

# Define groups
opts$groupA <- c(
  "E8.5_Dnmt3aKO_Dnmt3bWT", 
  "E8.5_Dnmt3aHET_Dnmt3bKO", 
  "E8.5_Dnmt3aHET_Dnmt3bWT", 
  "E8.5_Dnmt3aKO_Dnmt3bHET", 
  "E8.5_Dnmt3aKO_Dnmt3bKO", 
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

# Remove some hits
dt <- dt[gene!="Xist"]

##########
## Plot ##
##########

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
    # axis.text.x = element_blank(),
    axis.text.x = element_text(angle= -76 - 360 / length(unique(to.plot$celltype)) * seq_along(to.plot$celltype))
  )

# for (i in unique(dt$celltype)) {
#   to.plot <- dt[celltype==i] %>% .[!is.na(sig)] 
#   p <- gg_volcano_plot(to.plot, top_genes=15)
#   
#   pdf(sprintf("%s/%s_vs_%s_%s_volcano.pdf",io$outdir,opts$groupA,opts$groupB,i), width=9, height=5, useDingbats = F)
#   print(p)
#   dev.off()
# }



# to.plot <- dt %>% copy %>%
#   .[!is.na(logFC) & sig==T] %>%
#   .[,sign:=c("up in Mutant", "up in WT")[as.numeric(logFC>0)+1]] %>%
#   .[,.N,by=c("sign","celltype","groupA","groupB")]