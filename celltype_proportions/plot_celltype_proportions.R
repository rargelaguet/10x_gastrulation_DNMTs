here::i_am("celltype_proportions/plot_celltype_proportions.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
# p$add_argument('--samples',         type="character",       nargs="+",   help='Samples')
p$add_argument('--celltype_label', type="character", help='Cell type label')
p$add_argument('--outdir',          type="character",                               help='Output file')

args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
args$metadata <- file.path(io$basedir,"results/mapping/sample_metadata_after_mapping.txt.gz")
args$celltype_label <- "celltype.mapped"
args$outdir <- file.path(io$basedir,"results/celltype_proportions/test")
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings = F)
dir.create(file.path(args$outdir,"per_sample"), showWarnings = F)
dir.create(file.path(args$outdir,"per_class"), showWarnings = F)

# Options
opts$remove_ExE_cells <- TRUE

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  # .[pass_rnaQC==TRUE & doublet_call==FALSE & !is.na(eval(as.name(args$celltype_label)))] %>%
  .[pass_rnaQC==TRUE & !is.na(eval(as.name(args$celltype_label)))] %>%
  setnames(args$celltype_label,"celltype")

# Rename "_" to " " in cell types
# to.plot[,celltype:=stringr::str_replace_all(celltype,"_"," ")]
# names(opts$celltype.colors) <- names(opts$celltype.colors) %>% stringr::str_replace_all("_"," ")

if (opts$remove_ExE_cells) {
  print("Removing ExE cells...")
  sample_metadata <- sample_metadata %>%
    .[!celltype%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

################
## Per sample ##
################

to.plot <- sample_metadata %>%
  .[,N:=.N,by="alias"] %>%
  # .[,celltype:=stringr::str_replace_all(celltype,opts$aggregate.celltypes)] %>%
  .[,.(N=.N, celltype_proportion=.N/unique(N)),by=c("alias","celltype")] %>%
  setorder(alias)  %>% .[,alias:=factor(alias,levels=unname(opts$sample2alias))]

# Define colours and cell type order
opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors) %in% unique(to.plot$celltype)]
to.plot[,celltype:=factor(celltype, levels=rev(names(opts$celltype.colors)))]

samples.to.plot <- unique(to.plot$alias)

for (i in samples.to.plot) {
  p <- ggplot(to.plot[alias==i], aes(x=celltype, y=N)) +
    geom_bar(aes(fill=celltype), stat="identity", color="black") +
    scale_fill_manual(values=opts$celltype.colors) +
    facet_wrap(~alias, nrow=1, scales="fixed") +
    coord_flip() +
    labs(y="Number of cells") +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=rel(0.9)),
      axis.title.x = element_text(color="black", size=rel(0.9)),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1), color="black"),
      axis.text.x = element_text(size=rel(1), color="black")
    )
  
  pdf(file.path(args$outdir,sprintf("per_sample/celltype_proportions_%s_horizontal_barplots.pdf",i)), width=7, height=5)
  print(p)
  dev.off()
}

###############
## Per class ##
###############

to.plot <- sample_metadata %>%
  .[,N:=.N,by="class"] %>%
  .[,.(N=.N, celltype_proportion=.N/unique(N)),by=c("class","alias","celltype")] 

# Define colours and cell type order
opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors) %in% unique(to.plot$celltype)]
to.plot[,celltype:=factor(celltype, levels=rev(names(opts$celltype.colors)))]

classes.to.plot <- unique(to.plot$class)

for (i in classes.to.plot) {
  
  p <- ggplot(to.plot[class==i], aes(x=celltype, y=N)) +
    geom_bar(aes(fill=celltype), stat="identity", color="black") +
    scale_fill_manual(values=opts$celltype.colors) +
    facet_wrap(~alias, nrow=1, scales="fixed") +
    coord_flip() +
    labs(y="Number of cells") +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=rel(0.65)),
      axis.title.x = element_text(color="black", size=rel(0.9)),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(0.9), color="black"),
      axis.text.x = element_text(size=rel(0.75), color="black")
    )
  
  if (length(unique(to.plot[class==i,alias]))>=6) { width <- 15 } else { width <- 8 }
  
  pdf(file.path(args$outdir,sprintf("per_class/celltype_proportions_%s_horizontal_barplots.pdf",i)), width=width, height=4)
  print(p)
  dev.off()
}

######################
## Stacked barplots ##
######################

to.plot <- sample_metadata %>%
  .[,N:=.N,by="alias"] %>%
  .[,.(N=.N, celltype_proportion=.N/unique(N)),by=c("alias","celltype","class")] %>%
  setorder(alias)  %>% .[,alias:=factor(alias,levels=unname(opts$sample2alias))]

p <- ggplot(to.plot, aes(x=alias, y=celltype_proportion)) +
  geom_bar(aes(fill=celltype), stat="identity", color="black") +
  facet_wrap(~class, scales = "free_x", nrow=1) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    # axis.text.x = element_text(color="black", size=rel(0.75)),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

pdf(sprintf("%s/celltype_proportions_stacked_barplots.pdf",args$outdir), width=12, height=5)
print(p)
dev.off()

# Completion token
file.create(file.path(args$outdir,"completed.txt"))
