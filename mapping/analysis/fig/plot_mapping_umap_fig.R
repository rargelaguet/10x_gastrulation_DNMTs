here::i_am("mapping/analysis/plot_mapping_umap.R")

source(here::here("settings.R"))
source(here::here("mapping/analysis/plot_utils.R"))

io$query_metadata <- file.path(io$basedir,"results/mapping/sample_metadata_after_mapping.txt.gz")
io$atlas_metadata <- file.path(io$atlas.basedir,"sample_metadata.txt.gz")
io$outdir <- file.path(io$basedir,"results/mapping/pdf/fig"); dir.create(io$outdir, showWarnings = F)

#####################
## Define settings ##
#####################

# Options
opts$remove_ExE_cells <- FALSE
opts$subset_atlas <- TRUE

# Dot size
opts$size.mapped <- 0.18
opts$size.nomapped <- 0.1

# Transparency
opts$alpha.mapped <- 0.65
opts$alpha.nomapped <- 0.35

opts$samples <- c(
  "WT_1", "WT_2", "WT_3", "WT_4", "WT_6",
  "Dnmt3a_KO_1", "Dnmt3a_KO_12", "Dnmt3a_KO_2", "Dnmt3a_KO_Dnmt3b_HET_1", "Dnmt3a_KO_Dnmt3b_HET_2",
  "Dnmt3b_KO_1", "Dnmt3b_KO_2", "Dnmt3b_KO_6", "Dnmt3b_KO_7", "Dnmt3b_KO_9", 
  "Dnmt1_KO_10", "Dnmt1_KO_15", "Dnmt1_KO_2", "Dnmt1_KO_3", "Dnmt1_KO_9"
  # "Dnmt3ab_KO_1", "Dnmt3ab_KO_2"
)


#########################
## Load query metadata ##
#########################

sample_metadata <- fread(io$query_metadata) %>%
  .[,alias:=factor(alias,levels=opts$samples)] %>%
  .[class=="Dnmt3a_KO_Dnmt3b_HET",class:="Dnmt3a_KO"] %>%
  .[pass_rnaQC==TRUE & alias%in%opts$samples & !is.na(closest.cell)]

stopifnot("closest.cell"%in%colnames(sample_metadata))

if (opts$remove_ExE_cells) {
  print("Removing ExE cells...")
  sample_metadata <- sample_metadata %>% .[!celltype.mapped%in%c("ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

################
## Load atlas ##
################

# Load atlas cell metadata
meta_atlas <- fread(io$atlas_metadata) %>%
  # .[celltype%in%opts$celltypes] %>%
  .[stripped==F & doublet==F]

if (opts$remove_ExE_cells) {
  print("Removing ExE cells...")
  meta_atlas <- meta_atlas %>% .[!celltype%in%c("ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

if (opts$subset_atlas) {
  meta_atlas <- meta_atlas[sample.int(50000)]
}

# Extract precomputed dimensionality reduction coordinates
umap.dt <- meta_atlas %>%
  .[,c("cell","umapX","umapY","celltype")] %>%
  setnames(c("umapX","umapY"),c("V1","V2"))

##############################
## Define plotting function ##
##############################


###########################
## Plot multiple samples ##
###########################

ncells.to.subset <- 2500

to.plot <- opts$samples %>% map(function(i) {
  sample_metadata.subset <- sample_metadata[alias==i] %>% .[sample.int(min(ncells.to.subset,nrow(.)))]
  umap.dt %>% copy %>%
    .[,index:=match(cell, sample_metadata.subset[alias==i,closest.cell] )] %>% 
    .[,mapped:=as.factor(!is.na(index))] %>% 
    .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas","Query"))] %>%
    .[,sample:=factor(i,levels=opts$samples)] %>%
    setorder(mapped) 
}) %>% rbindlist

p <- plot.dimred(to.plot, query.label = "Query", atlas.label = "Atlas") + 
  facet_wrap(~sample, nrow=4) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(color="black", size=rel(0.8))
  )
  
pdf(file.path(io$outdir,"umap_mapped_all_samples.pdf"), width=13, height=13)
print(p)
dev.off()

###############################
## Plot one class at a time ##
###############################

classes.to.plot <- unique(sample_metadata$class)

# Subset classes to similar number of cells
# ncells.to.subset <- min(sample_metadata[,.N,by="class"][["N"]])
ncells.to.subset <- 5000

for (i in classes.to.plot) {
  
  sample_metadata.subset <- sample_metadata[class==i] %>% .[sample.int(min(ncells.to.subset,nrow(.)))]
  to.plot <- umap.dt %>% copy %>%
    .[,index:=match(cell, sample_metadata.subset[,closest.cell] )] %>% 
    .[,mapped:=as.factor(!is.na(index))] %>% 
    .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas",i))] %>%
    setorder(mapped) 
  
  p <- plot.dimred(to.plot, query.label = i, atlas.label = "Atlas") + theme(legend.position = "none")
  
  pdf(sprintf("%s/per_class/umap_mapped_%s.pdf",io$outdir,i), width=8, height=6.5)
  print(p)
  dev.off()
}

#############################
## Plot WT and KO together ##
#############################

# Subsample query cells to have the same N per class
# sample_metadata_subset <- sample_metadata %>% .[,.SD[sample.int(n=.N, size=4500)], by=c("stage","class2")]

# i <- "E7.5"
for (i in args$classes) {
  
  to.plot <- umap.dt %>% copy %>%
    .[,index.wt:=match(cell, sample_metadata[class=="WT" & stage==i,closest.cell] )] %>%
    .[,index.ko:=match(cell, sample_metadata[class==i & stage==i,closest.cell] )] %>%
    .[,mapped.wt:=c(0,-10)[as.numeric(as.factor(!is.na(index.wt)))]] %>%
    .[,mapped.ko:=c(0,10)[as.numeric(as.factor(!is.na(index.ko)))]] %>%
    .[,mapped:=factor(mapped.wt + mapped.ko, levels=c("0","-10","10"))] %>%
    .[,mapped:=plyr::mapvalues(mapped, from = c("0","-10","10"), to = c("Atlas","WT",i))] %>% setorder(mapped)
  
  p <- plot.dimred.wtko(to.plot, wt.label = "WT", ko.label = i, nomapped.label = "Atlas") +
    theme(legend.position = "top", axis.line = element_blank())
  
  pdf(sprintf("%s/per_class/umap_mapped_%s_WT_and_KO.pdf",io$outdir,i), width=5.5, height=6.5)
  print(p)
  dev.off()
}
# Completion token
file.create(file.path(io$outdir,"completed.txt"))
