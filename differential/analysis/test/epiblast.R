##############
## Settings ##
##############

# Load settings
source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
source("/Users/ricard/10x_gastrulation_DNMTs/differential/analysis/utils.R")

# I/O
io$diff.dir <- paste0(io$basedir,"/results/differential")
io$outdir <- paste0(io$basedir,"/results/differential/test")

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

opts$celltypes = c(
  "Epiblast"
)

###############################
## Load Nueral crest markers ##
###############################

marker_genes.dt <- fread(io$atlas.marker_genes)
# neural_crest.dt <- fread("/Users/ricard/data/gastrulation10x/results/differential/celltypes/E8.5/Neural_crest_vs_Forebrain_Midbrain_Hindbrain.txt.gz") %>%
#   .[sig==T & logFC<0]

#############################################
## Load results from differential analysis ##
#############################################

diff.dt <- opts$groupA %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- sprintf("%s/%s_vs_%s_%s.txt.gz", io$diff.dir,i,opts$groupB,j)
  if (file.exists(file)) fread(file) %>% .[,c(1,2,4,6,7,10,11)] %>% 
    # setnames(c(sprintf("N_%s",i),sprintf("N_%s",opts$groupB)),c("N_groupA","N_groupB")) %>%
    .[,c("celltype","groupA","groupB"):=list(j,i,opts$groupB)]
}) %>% rbindlist }) %>% rbindlist

unique(diff.dt$groupA)
####################
## Filter results ##
####################

# select neural crest markers
diff_filt.dt <- diff.dt[gene%in%unique(marker_genes.dt$gene)]

##########
## Plot ##
##########

for (i in unique(diff_filt.dt$groupA)) {
  to.plot <- diff_filt.dt[groupA==i]
  
  negative_hits <- to.plot[sig==TRUE & logFC<0,gene]
  positive_hits <- to.plot[sig==TRUE & logFC>0,gene]
  all <- nrow(to.plot)
  
  xlim <- max(abs(to.plot$logFC), na.rm=T)
  ylim <- max(-log10(to.plot$padj_fdr+1e-100), na.rm=T)
  
  p <- ggplot(to.plot, aes(x=logFC, y=-log10(padj_fdr+1e-100))) +
    labs(x="Log fold change", y=expression(paste("-log"[10],"(p.value)"))) +
    geom_segment(aes(x=0, xend=0, y=0, yend=ylim-1), color="orange", size=0.5) +
    ggrastr::geom_point_rast(aes(color=sig, size=sig)) +
    scale_color_manual(values=c("black","red")) +
    scale_size_manual(values=c(0.75,1.25)) +
    scale_x_continuous(limits=c(-xlim-1.5,xlim+1.5)) +
    scale_y_continuous(limits=c(0,ylim+3)) +
    annotate("text", x=-xlim-0.5, y=0, size=3, label=sprintf("Up in %s (N=%s)",i,length(negative_hits))) +
    annotate("text", x=xlim+0.5, y=0, size=3, label=sprintf("Up in %s (N=%s)",opts$groupB,length(positive_hits))) +
    ggrepel::geom_text_repel(data=to.plot[sig==T], aes(x=logFC, y=-log10(padj_fdr+1e-100), label=gene), size=4, max.overlaps=100) +
    theme_classic() +
    theme(
      axis.text = element_text(size=rel(0.75), color='black'),
      axis.title = element_text(size=rel(1.0), color='black'),
      legend.position="none"
    )
  
  
  pdf(sprintf("%s/%s_vs_%s_volcano.pdf",io$outdir,i,opts$groupB), width=9, height=5)
  print(p)
  dev.off()
}


fwrite(diff_filt.dt[!is.na(logFC)], file="/Users/ricard/data/10x_gastrulation_DNMTs/results/differential/test/diff_epiblast_stephen.txt.gz", sep="\t", na="NA")
