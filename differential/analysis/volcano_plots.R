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
io$indir <- file.path(io$basedir,"results/differential")
io$outdir <- file.path(io$basedir,"results/differential/pdf/volcano_plots"); dir.create(io$outdir, showWarnings = F, recursive = T)

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

# Remove some hits
# diff.dt <- diff.dt[gene!="Xist"]
# diff.dt <- diff.dt[!grepl("mt-",gene)]
# diff.dt <- diff.dt[!grepl("Rps|Rpl",gene)]
# diff.dt <- diff.dt[!grepl("Rik",gene)]
# diff.dt <- diff.dt[!grepl("^Hb",gene)]
# diff.dt <- diff.dt[!gene%in%fread(io$gene_metadata)[chr=="chrY",symbol]]

# Filter by minimum number of cells per group
# opts$min.cells <- 30
# diff.dt <- diff.dt[groupA_N>opts$min.cells & groupB_N>opts$min.cells]

# Remove hits that are differentially expressed in all cell type comparisons
# foo <- diff.dt[,mean(sig),by=c("gene")] %>% .[V1>0] %>% setorder(-V1)

# Subset to lineage markers
marker_genes.dt <- fread(io$atlas.marker_genes)
diff_markers.dt <- diff.dt[gene%in%unique(marker_genes.dt$gene)]

##########
## Plot ##
##########

# i <- opts$ko.classes[1]; j <- opts$celltypes[1]
for (i in opts$ko.classes) {
  outdir <- file.path(io$outdir,i); dir.create(outdir, showWarnings = F)
  celltypes.to.plot <- unique(diff_markers.dt[class==i,celltype]) %>% as.character
  for (j in celltypes.to.plot) {
    to.plot <- diff_markers.dt[celltype==j & class==i] %>% .[!is.na(sig)] 
    if (nrow(to.plot)>0) {
      p <- gg_volcano_plot(to.plot, top_genes = 35, groupA = opts$wt.class, groupB = i)
      
      pdf(file.path(outdir,sprintf("%s_vs_%s_%s_volcano.pdf",i,opts$wt.class,j)), width=9, height=6)
      print(p)
      dev.off()
    }
  }
}
