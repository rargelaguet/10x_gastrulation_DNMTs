#####################
## Define settings ##
#####################

# load default setings
source(here::here("settings.R"))
source(here::here("utils.R"))

# Options
opts$TFs <- fread(file.path(io$basedir,"shiny/TFs.txt"), header = F)[[1]] %>% str_to_title

# I/O
io$rna.diff.pseudobulk <- file.path(io$basedir,"results/rna/differential/pseudobulk")
io$outdir <- file.path(io$basedir,"results/rna/differential/pseudobulk/pdf")

##############################
## Load precomputed results ##
##############################

# source(here::here("rna_atac/load_rna_atac_pseudobulk.R"))

diff.dt <- (1:length(opts$celltypes)) %>% map(function(i) {
  (i:length(opts$celltypes)) %>% map(function(j) {
    if (i!=j) {
      fread(sprintf("%s/%s_vs_%s.txt.gz", io$rna.diff.pseudobulk,opts$celltypes[[i]],opts$celltypes[[j]])) %>%
        .[gene%in%opts$TFs] %>% .[,gene:=toupper(gene)]  %>%
        .[,diff:=round(groupB-groupA,2)] %>% .[,c("groupA","groupB"):=NULL] %>%
        .[,c("groupA","groupB"):=list(opts$celltypes[[i]],opts$celltypes[[j]])] %>%
        sort.abs("diff") 
    }
  }) %>% rbindlist
}) %>% rbindlist

# save
fwrite(diff.dt, file.path(io$outdir,"diff_rna_pseudobulk.txt.gz"))

##########
## Plot ##
##########

i <- "Surface_ectoderm"
j <- "Erythroid3"
# celltypes.to.plot <- c("Gut","Erythroid3")
# genes.to.plot <- c("TAL1")

to.plot <- diff.dt[groupA==i & groupB==j] %>% .[,gene:=factor(gene,levels=rev(gene))]
to.plot <- diff.dt[groupB==i & groupA==j] %>% .[,gene:=factor(gene,levels=rev(gene))]

p <- ggplot(to.plot, aes(x=gene, y=diff)) +
  geom_point(aes(color=abs(diff), alpha=abs(diff))) +
  ggrepel::geom_text_repel(data=head(to.plot[diff>0],n=10), aes(x=gene, y=diff, label=gene), size=10) +
  scale_color_gradient(low = "gray80", high = "red") +
  scale_alpha_continuous(range=c(0.25,1)) +
  theme_classic() +
  labs(y="Differential RNA expression", x="") +
  annotate("text", x=35, y=-12, size=4, label=sprintf("(+) %s",i)) +
  annotate("text", x=35, y=12, size=4, label=sprintf("(+) %s",j)) +
  coord_flip() +
  geom_segment(x=0, xend=800, y=0, yend=0, color="black", size=0.25, linetype="dashed") +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size=rel(1.0), color="black")
  )
