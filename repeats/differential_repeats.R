source(here::here("settings.R"))
source(here::here("utils.R"))


##############
## Settings ##
##############

# Load settings
source(here::here("differential/analysis/utils.R"))

# I/O
io$indir <- paste0(io$basedir,"/results/differential")
io$outdir <- paste0(io$basedir,"/results/differential/repeats"); dir.create(io$outdir, showWarnings = F)

# Options

# Define classes to plot
opts$classes <- c(
  "WT", 
  "Dnmt3a_KO", 
  "Dnmt3b_KO",
  "Dnmt1_KO"
  # "Dnmt3ab_KO"
)

# Rename cell types
opts$rename_celltypes <- c(
  "Erythroid3" = "Erythroid",
  "Erythroid2" = "Erythroid",
  "Erythroid1" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors"
)

stringr::str_replace_all(opts$celltypes, opts$rename_celltypes) %>% unique()

opts$celltypes <- c(
  "Caudal_epiblast",
  "Notochord", 
  "Def._endoderm",
  "Gut",
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
  "Blood_progenitors",
  "Erythroid",
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
)

opts$repeat_classes <- c(
  "LINE_L1",
  "LINE_L2",
  "LTR_ERV1",
  "LTR_ERVK",
  "LTR_ERVL",
  "LTR_MaLR",
  "major_satellite",
  "minor_satellite",
  "rRNA",           
  "SINE_Alu_B1",
  "SINE_B2",
  "SINE_B4",
  "IAP"  
)

opts$min.exp <- 1 # exclude celltypes or repeat classes with universally low expression (i.e. in both WT and KO)

opts$min.cells <- 50


# load metadata for filtering out celltypes with low numbers of cells

meta <- fread(io$metadata) %>% 
  .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,opts$rename_celltypes)] %>% 
  .[, barcode2 := gsub("\\-.*", "", barcode)]

n_celltypes <- meta[!is.na(celltype.mapped) & dataset == "KO", .N, .(celltype.mapped, class)] %>% 
  dcast(celltype.mapped ~ class, value.var = "N") %>% 
  .[, Dnmt1 := ifelse(WT < Dnmt1_KO, WT, Dnmt1_KO)] %>% 
  .[, Dnmt3a := ifelse(WT < Dnmt3a_KO, WT, Dnmt3a_KO)] %>% 
  .[, Dnmt3b := ifelse(WT < Dnmt3b_KO, WT, Dnmt3b_KO)] %>% 
  melt(id.vars = "celltype.mapped", measure.vars = c("Dnmt1", "Dnmt3a", "Dnmt3b"), variable.name = "class", value.name = "min_cells")

######################
## Load repeat data ##
######################

exp <- fread(paste0(io$basedir, "/repeats/repeats_expr.txt.gz")) %>% 
  .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,opts$rename_celltypes)] %>% 
  .[celltype.mapped %in% celltypes] %>% 
  .[repeat_class %in% opts$repeat_classes] %>% 
  setkey(alias, celltype.mapped)

# code to double check expression values are correct
# exp_raw <- dir("~/data/10x_gastrulation_DNMTs/repeats/10X_cellwise_repeat_content/", full = TRUE) %>% 
#   set_names(basename(.) %>%strsplit("_") %>% map_chr(~paste(.x[1], .x[2], .x[3], sep = "_")) %>% gsub("_CB", "", .) ) %>% 
#   map2(names(.), ~fread(.x)[, alias := .y]) %>% 
#   rbindlist() %>% 
#   setnames(colnames(.)[1 : (ncol(.) - 1)], c(colnames(.)[2 : (ncol(.) - 1)], "remove")) %>% 
#   .[, barcode2 := gsub("CB:", "", `cell barcode`)] %>% 
#   merge(meta[, .(alias, barcode2, class, celltype.mapped)], by = c("alias", "barcode2"))
# 
# exp_by_alias <- exp_raw[, map(.SD, sum), by = c("alias", "celltype.mapped"), .SDcol = c(opts$repeat_classes, "total reads")] %>% 
#   melt(id.vars = c("alias", "celltype.mapped", "total reads")) %>% 
#   .[!is.na(celltype.mapped)] %>% 
#   setkey(alias, celltype.mapped) %>% 
#   .[, expr := log2(1e6 * value / `total reads`)]

# exp_raw_long <- melt(exp_raw[, .SD, .SDcol = c("barcode2","alias", "celltype.mapped", "total reads", opts$repeat_classes)], id.vars = c("barcode2","alias", "celltype.mapped", "total reads")) %>% 
#     .[!is.na(celltype.mapped)] %>%
#     setkey(alias, celltype.mapped) %>%
#     .[, expr := log2(1e6 * value / `total reads`)]
# 
# exp_raw_long[alias == "Dnmt1_KO_1" & celltype.mapped == "NMP", summary(expr) %>% as.list, variable]
# exp_raw_long[alias == "WT_2" & celltype.mapped == "NMP", summary(expr) %>% as.list, variable]
# exp_raw_long[alias == "Dnmt1_KO_1" & celltype.mapped == "ExE_ectoderm", summary(expr) %>% as.list, variable]
# exp_raw_long[alias == "WT_2" & celltype.mapped == "ExE_ectoderm", summary(expr) %>% as.list, variable]
# 
# 
# exp_per_class[class == "Dnmt1_KO" & celltype.mapped == "NMP"]
# exp_per_class[class == "WT" & celltype.mapped == "NMP"]
# exp_per_class[class == "Dnmt1_KO" & celltype.mapped == "NMP"]
# exp_per_class[class == "Dnmt1_KO" & celltype.mapped == "NMP"]


# compute expression after merging by class (i.e. WT, Dnmt1_KO etc.)

exp_per_class <- copy(exp) %>% 
  .[, .(nreads = sum(nreads), total_reads = sum(total_reads)), .(celltype.mapped, class,  repeat_class) ] %>% 
  .[, expr := log2(1e6 * nreads/total_reads)]

diff.dt <- dcast(exp_per_class, celltype.mapped + repeat_class ~ class, value.var = "expr") %>% 
  .[, Dnmt1 := Dnmt1_KO - WT] %>% 
  .[, Dnmt3a := Dnmt3a_KO - WT] %>% 
  .[, Dnmt3b := Dnmt3b_KO - WT] %>% 
  .[, D1_max := ifelse(Dnmt1_KO > WT, Dnmt1_KO, WT)] %>% 
  .[, D3a_max := ifelse(Dnmt3a_KO > WT, Dnmt3a_KO, WT)] %>% 
  .[, D3b_max := ifelse(Dnmt3b_KO > WT, Dnmt3b_KO, WT)] %>% 
  melt(
    id.vars = c("celltype.mapped", "repeat_class"), 
    measure.vars = list(c("Dnmt1", "Dnmt3a", "Dnmt3b"), c("D1_max", "D3a_max", "D3b_max")),
    value.name = c("logFC", "max_expr"),
    variable.name = "class"
    ) %>% 
  .[class == 1, class := "Dnmt1"] %>% 
  .[class == 2, class := "Dnmt3a"] %>% 
  .[class == 3, class := "Dnmt3b"] %>% 
  .[, celltype := celltype.mapped] %>% 
  .[, gene := repeat_class]


# filter for number of cells

diff.dt.filt <- merge(diff.dt, n_celltypes, by = c("celltype.mapped", "class")) %>% 
  .[min_cells >= opts$min.cells]

# optionally, we might want to keep all plots the same with respect to which celltypes are plotted
celltype_filter <- n_celltypes[min_cells >= opts$min.cells]  %>% 
  .[, .N, celltype.mapped] %>% 
  .[N == 3, celltype.mapped]

opts$celltypes <- opts$celltypes[opts$celltypes %in% celltype_filter]

# # filter out lowly expressed
diff.dt.filt[max_expr < opts$min.exp, logFC := 0]

###################################
## Heatmaps, one class at a time ##
###################################



opts$max.logFC <- 15; opts$min.logFC <- -15
opts$ko.classes <- rev(diff.dt[,unique(class)])



for (i in opts$ko.classes) {
  
  tmp <- diff.dt.filt[class==i] %>%
    .[logFC>opts$max.logFC,logFC:=opts$max.logFC] %>%
    .[logFC<opts$min.logFC,logFC:=opts$min.logFC]
  
  to.plot <- expand.grid(X = unique(tmp$gene), Y = unique(tmp$celltype)) %>% 
    as.data.table %>% setnames(c("gene","celltype")) %>%
    merge(tmp, by=c("gene","celltype"), all.x=T)
  
  # # filter out celltypes with low expression of all repeats
  # to.plot[, allna := all(is.na(logFC)), celltype]
  # (remove <- to.plot[allna==T, unique(celltype.mapped)])
  # to.plot <- to.plot[!celltype.mapped %in% remove]
  # 
  # # filter out repeat classes with low expression in all celltypes
  # to.plot[, allna := all(is.na(logFC)), gene]
  # (remove <- to.plot[allna==T, unique(gene)])
  # to.plot <- to.plot[!gene %in% remove]
  
  # reorder celltyeps
  to.plot <- to.plot[, celltype := factor(celltype, levels = rev(opts$celltypes))] %>% 
    .[!is.na(celltype)]
  
  # format repeat name
  to.plot[, gene := gsub("_", "\n", gene)]
  
  p <- ggplot(to.plot, aes(x=gene, y=celltype, fill=logFC)) +
    geom_tile(color="black") +
    # viridis::scale_fill_viridis() +
    # scale_fill_gradientn(colours = terrain.colors(10)) +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, limits=c(opts$min.logFC-0.01,opts$max.logFC+0.01), na.value = 'grey') +
    labs(x="", y="") +
    #guides(x = guide_axis(angle = 90)) +
    theme_classic() +
    theme(
      legend.position = "top",
      axis.line = element_blank(),
      axis.text = element_text(color="black", size=rel(0.75))
    )
  p
  pdf(sprintf("%s/%s_DE_heatmap_repeats.pdf",io$outdir,i), width=12, height=7)
  print(p)
  dev.off()
}
p
