# here::i_am("processing/1_create_seurat_rna.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

# suppressPackageStartupMessages(library(Seurat))

#####################
## Define settings ##
#####################

io$repeats_expr <- file.path(io$basedir,"results/repeats/repeats_expr.txt.gz")
io$outdir <- file.path(io$basedir,"results/repeats/pdf"); dir.create(io$outdir, showWarnings = F)

opts$samples <- c(
  "Dnmt1_KO_1",
  "Dnmt1_KO_2",
  "Dnmt1_KO_3",
  "Dnmt3a_KO_13",
  "Dnmt3a_KO_14",
  "Dnmt3a_KO_1",
  "Dnmt3a_KO_2",
  "Dnmt3b_KO_11",
  "Dnmt3b_KO_12",
  "Dnmt3b_KO_1",
  "Dnmt3b_KO_2",
  "WT_1",
  "WT_2",
  "WT_3",
  "WT_4",
  "WT_5",
  "WT_6",
  "WT_7"
)

########################
## Load cell metadata ##
########################

cell_metadata.dt <- fread(io$metadata) %>% 
  .[pass_rnaQC==TRUE & alias%in%opts$samples]
table(cell_metadata.dt$alias)

class_metadata.dt <- cell_metadata.dt[,.N,by=c("class","celltype.mapped")] %>% .[N>=100]

########################
## Load TE expression ##
########################

repeats.dt <- fread(io$repeats_expr)

table(repeats.dt$repeat_class)

##########
## Plot ##
##########

celltypes.to.plot <- c("ExE_endoderm","Mesenchyme","Haematoendothelial_progenitors","Surface_ectoderm")
celltypes.to.plot <- opts$celltypes
elements.to.plot <- unique(repeats.dt$repeat_class)

for (i in elements.to.plot) {
  
  to.plot <- repeats.dt[repeat_class==i & celltype.mapped%in%celltypes.to.plot] %>% 
    merge(class_metadata.dt,by=c("class","celltype.mapped"))
  
  give_n <- function(x){ return(c(y = max(to.plot$expr), label = length(x))) }
  
  my_comparisons <- list( c("Dnmt1_KO", "WT"))
  
  p <- ggplot(to.plot, aes(x=class, y=expr)) +
    geom_boxplot(aes(fill = class), alpha=0.75) +
    geom_jitter(aes(fill=class), shape=21, size=1.5, alpha=0.75, width=0.05, stroke=0.15, show.legend = F) +
    # facet_wrap(~celltype.mapped, scales="fixed", nrow=1) +
    facet_wrap(~celltype.mapped, scales="fixed") +
    scale_fill_manual(values=opts$classes.colors[unique(to.plot$class)]) +
    # stat_compare_means(aes(label = as_name(paste0("p = ", ..p.format..))), comparisons = my_comparisons, method="t.test", size=3) +
    theme_classic() +
    labs(x="", y="RNA expression (log2)") +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      axis.text.y = element_text(color="black", size=rel(0.8)),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  # pdf(file.path(io$outdir,sprintf("%s_repeats_boxplots.pdf",i)), width=8, height=3)
  pdf(file.path(io$outdir,sprintf("%s_repeats_boxplots.pdf",i)), width=10, height=10)
  print(p)
  dev.off()
}

