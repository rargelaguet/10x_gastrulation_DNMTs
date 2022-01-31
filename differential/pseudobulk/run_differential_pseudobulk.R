source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$sce.pseudobulk <- file.path(io$basedir,"results/pseudobulk/SingleCellExperiment_pseudobulk_class_celltype.rds")
io$outdir <- file.path(io$basedir,"results/differential/pseudobulk"); dir.create(io$outdir, showWarnings = F)

# Options
opts$min.cells <- 25

opts$ko.classes <- c(
  "Dnmt3a_KO", 
  # "Dnmt3a_HET_Dnmt3b_KO", 
  # "Dnmt3a_HET_Dnmt3b_WT", 
  # "Dnmt3a_KO_Dnmt3b_HET", 
  "Dnmt3ab_KO", 
  "Dnmt3b_KO",
  "Dnmt1_KO"
)

opts$wt.class <- "WT"

##############################
## Load pseudobulk RNA data ##
##############################

sce <- readRDS(io$sce.pseudobulk)

# Add metadata
sce$class <- stringr::str_split(colnames(sce), pattern = "-") %>% map_chr(1)
sce$celltype <- stringr::str_split(colnames(sce), pattern = "-") %>% map_chr(2)

# Filter by minimum number of cells (this has to be done first)
sce <- sce[,names(which(metadata(sce)$n_cells>=opts$min.cells))]

# Subset classes
sce <- sce[,sce$class%in%opts$classes]
metadata(sce)$n_cells <- metadata(sce)$n_cells[stringr::str_split(names(metadata(sce)$n_cells), pattern = "-") %>% map_chr(1) %in% c(opts$ko.classes,opts$wt.class)]

# Subset celltypes
# sce <- sce[,sce$celltype%in%opts$celltypes]
# metadata(sce)$n_cells <- metadata(sce)$n_cells[stringr::str_split(names(metadata(sce)$n_cells), pattern = "-") %>% map_chr(2) %in% opts$celltypes]

####################################################
## Differential expression per class and celltype ##
####################################################

# Subset WT samples from the SingleCellExperiment
sce.wt <- sce[,sce$class==opts$wt.class]
metadata(sce.wt)$n_cells <- metadata(sce.wt)$n_cells[stringr::str_split(names(metadata(sce.wt)$n_cells), pattern = "-") %>% map_chr(1) %in% opts$wt.class]
names(metadata(sce.wt)$n_cells) <- stringr::str_split(names(metadata(sce.wt)$n_cells), pattern = "-") %>% map_chr(2)
colnames(sce.wt) <- sce.wt$celltype

# Select celltypes with sufficient number of cells
celltypes.wt <- names(which(metadata(sce.wt)$n_cells >= opts$min.cells))

# i <- "Dnmt1_KO"
for (i in opts$ko.classes) {
  
  # Subset KO samples from the SingleCellExperiment
  sce.ko <- sce[,sce$class==i]
  metadata(sce.ko)$n_cells <- metadata(sce.ko)$n_cells[stringr::str_split(names(metadata(sce.ko)$n_cells), pattern = "-") %>% map_chr(1) == i]
  names(metadata(sce.ko)$n_cells) <- stringr::str_split(names(metadata(sce.ko)$n_cells), pattern = "-") %>% map_chr(2)
  colnames(sce.ko) <- sce.ko$celltype
  
  # Select celltypes with sufficient number of cells
  celltypes.ko <- names(which(metadata(sce.ko)$n_cells >= opts$min.cells))
  
  celltypes.to.use <- intersect(celltypes.wt,celltypes.ko)
  for (j in celltypes.to.use) {
      
    tmp <- data.table(
      gene = rownames(sce.ko),
      expr_ko = logcounts(sce.ko[,j])[,1] %>% round(2),
      expr_wt = logcounts(sce.wt[,j])[,1] %>% round(2)
    ) %>% .[,diff:=round(expr_ko-expr_wt,2)] %>% sort.abs("diff") 

    # save      
    outfile <- sprintf("%s/%s/%s_%s_vs_%s.txt.gz", io$outdir,i,j,opts$wt.class,i); dir.create(dirname(outfile), showWarnings = F)
    fwrite(tmp, outfile, sep="\t")
  }
}

#######################################
## Differential expression per class ##
#######################################

# # Subset WT samples from the SingleCellExperiment
# sce.wt <- sce[,sce$class==opts$wt.class]
# 
# for (i in opts$ko.classes) {
#   
#   # Subset KO samples from the SingleCellExperiment
#   sce.ko <- sce[,sce$class==i]
#   
#   tmp <- data.table(
#     gene = rownames(sce.ko),
#     expr_ko = logcounts(sce.ko)[,1] %>% round(2),
#     expr_wt = logcounts(sce.wt)[,1] %>% round(2)
#   ) %>% .[,diff:=round(expr_ko-expr_wt,2)] %>% sort.abs("diff") 
#     
#   # save      
#   outfile <- sprintf("%s/%s/%s_vs_%s.txt.gz", io$outdir,i,opts$wt.class,i); dir.create(dirname(outfile), showWarnings = F)
#   fwrite(tmp, outfile, sep="\t")
# }
