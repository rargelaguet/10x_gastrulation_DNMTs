###################################################################
## Plot dimensionality reduction of EB cells mapped to the atlas ##
###################################################################

source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
source("/Users/ricard/10x_gastrulation_DNMTs/utils.R")
source("/Users/ricard/10x_gastrulation_DNMTs/mapping/analysis/plot_utils.R")

io$outdir <- paste0(io$basedir,"/results/mapping/pdf/umap_overlayed_atlas")
dir.create(paste0(io$outdir,"/per_sample"),showWarnings = F)
dir.create(paste0(io$outdir,"/per_class"),showWarnings = F)

#############
## Options ##
#############

opts$classes <- c(
  "E8.5_Dnmt3aKO_Dnmt3bWT",
  "E8.5_WT",
  "E8.5_Dnmt3aHET_Dnmt3bKO",
  "E8.5_Dnmt3aHET_Dnmt3bWT",
  "E8.5_Dnmt3aKO_Dnmt3bHET",
  "E8.5_Dnmt3aKO_Dnmt3bKO",
  "E8.5_Dnmt3aWT_Dnmt3bKO",
  "E8.5_Dnmt1KO"
)

opts$wt.class <- "E8.5_WT"

# Dot size
opts$size.mapped <- 0.25
opts$size.nomapped <- 0.1

# Transparency
opts$alpha.mapped <- 0.8
opts$alpha.nomapped <- 0.30

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  .[pass_QC==TRUE & celltype.mapped%in%opts$celltypes & class%in%opts$classes] %>%
  .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped," ","_")] %>%
  .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,"/","_")]

# Rename samplees
foo <- sample_metadata[,c("batch","class")] %>% unique %>% .[,sample:=paste(class,1:.N,sep="_"),by="class"]
sample_metadata <- sample_metadata %>% merge(foo,by=c("batch","class"))

####################
## Load 10x atlas ##
####################

# Load atlas cell metadata
meta_atlas <- fread(io$atlas.metadata) %>%
  .[stripped==F & doublet==F]

# Extract precomputed dimensionality reduction coordinates
umap.dt <- meta_atlas[,c("cell","umapX","umapY","celltype")] %>%
  setnames(c("umapX","umapY"),c("V1","V2"))

#########################################################
## Plot dimensionality reduction: one embryo at a time ##
#########################################################

opts$samples <- unique(sample_metadata$sample)

for (i in opts$samples) {
  
  # Alone
  to.plot <- umap.dt %>% copy %>%
    .[,index:=match(cell, sample_metadata[sample==i,closest.cell] )] %>% 
    .[,mapped:=as.factor(!is.na(index))] %>% 
    .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas",i))] %>%
    setorder(mapped) 
  
  p <- plot.dimred(to.plot, query.label = i, atlas.label = "Atlas")
  
  pdf(sprintf("%s/per_sample/umap_mapped_%s.pdf",io$outdir,i), width=8, height=6.5)
  print(p)
  dev.off()
  
  # Together with WT
  if (sample_metadata[sample==i,class]!=opts$wt.class) {
    to.plot <- umap.dt %>% copy %>%
      .[,index.wt:=match(cell, sample_metadata[class==opts$wt.class,closest.cell] )] %>%
      .[,index.ko:=match(cell, sample_metadata[sample==i,closest.cell] )] %>%
      .[,mapped.wt:=c(0,-10)[as.numeric(as.factor(!is.na(index.wt)))]] %>%
      .[,mapped.ko:=c(0,10)[as.numeric(as.factor(!is.na(index.ko)))]] %>%
      .[,mapped:=factor(mapped.wt + mapped.ko, levels=c("0","-10","10"))] %>%
      .[,mapped:=plyr::mapvalues(mapped, from = c("0","-10","10"), to = c("Atlas",opts$wt.class,i))] %>% setorder(mapped)
    
    p <- plot.dimred.wtko(to.plot, wt.label = opts$wt.class, ko.label = i, nomapped.label = "Atlas")
    
    pdf(sprintf("%s/per_sample/umap_mapped_%s_withWT.pdf",io$outdir,i), width=8, height=6.5)
    print(p)
    dev.off()
  }
}

########################################################
## Plot dimensionality reduction: one class at a time ##
########################################################

for (i in opts$classes) {
  
  # Alone
  to.plot <- umap.dt %>% copy %>%
    .[,index:=match(cell, sample_metadata[class==i,closest.cell] )] %>% 
    .[,mapped:=as.factor(!is.na(index))] %>% 
    .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas",i))] %>%
    setorder(mapped) 
  
  p <- plot.dimred(to.plot, query.label = i, atlas.label = "Atlas")
  
  pdf(sprintf("%s/per_class/umap_mapped_%s.pdf",io$outdir,i), width=8, height=6.5)
  print(p)
  dev.off()
  
  
  # Together with WT
  if (i!=opts$wt.class) {
    to.plot <- umap.dt %>% copy %>%
      .[,index.wt:=match(cell, sample_metadata[class==opts$wt.class,closest.cell] )] %>%
      .[,index.ko:=match(cell, sample_metadata[class==i,closest.cell] )] %>%
      .[,mapped.wt:=c(0,-10)[as.numeric(as.factor(!is.na(index.wt)))]] %>%
      .[,mapped.ko:=c(0,10)[as.numeric(as.factor(!is.na(index.ko)))]] %>%
      .[,mapped:=factor(mapped.wt + mapped.ko, levels=c("0","-10","10"))] %>%
      .[,mapped:=plyr::mapvalues(mapped, from = c("0","-10","10"), to = c("Atlas",opts$wt.class,i))] %>% setorder(mapped)
    
    p <- plot.dimred.wtko(to.plot, wt.label = opts$wt.class, ko.label = i, nomapped.label = "Atlas")
    pdf(sprintf("%s/per_class/umap_mapped_%s_withWT.pdf",io$outdir,i), width=8, height=6.5)
    print(p)
    dev.off()
  }
}
