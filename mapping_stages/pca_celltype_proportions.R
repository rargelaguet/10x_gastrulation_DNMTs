
#####################
## Define settings ##
#####################

source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
source("/Users/ricard/10x_gastrulation_DNMTs/utils.R")

# I/O
io$outdir <- paste0(io$basedir,"/results/mapping_stages")

# Options
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

opts$remove.ExE.celltypes <- TRUE

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  .[pass_QC==TRUE & class%in%opts$classes] %>%
  .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped," ","_")] %>%
  .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,"/","_")]

if (opts$remove.ExE.celltypes) {
  sample_metadata <- sample_metadata %>%
    .[!celltype.mapped%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}
table(sample_metadata$class)

################
## Load query ##
################

dt.query <- sample_metadata %>% copy %>%
  .[,N:=.N,by="batch"] %>%
  setnames("celltype.mapped","celltype") %>%
  .[,.(celltype_proportion=.N/unique(N)),by=c("batch","celltype")]

sample_metadata.query <- sample_metadata %>%
  .[,.(nCount_RNA=mean(nCount_RNA)),c("batch","stage","class")] %>%
  .[,dataset:="Query"]

################
## Load atlas ##
################

sample_metadata.atlas <- fread(io$atlas.metadata) %>% 
  .[stage!="mixed_gastrulation"] %>%
  setnames("sample","batch") %>%
  .[,.(nCount_RNA=mean(nCount_RNA)),c("batch","stage")] %>%
  .[,c("class","dataset"):="Atlas"]

# Load cell type proportions in the atlas
dt.atlas <- fread(paste0(io$atlas.basedir,"/results/celltype_proportions/celltype_proportions.txt.gz")) %>%
  .[sample%in%sample_metadata.atlas$batch] %>%
  setnames("sample","batch")


#################
## Concatenate ##
#################

sample_metadata <- rbind(
  sample_metadata.query[,c("batch", "stage", "nCount_RNA", "class", "dataset")], 
  sample_metadata.atlas[,c("batch", "stage", "nCount_RNA", "class", "dataset")]
)

dt <- rbind(
  dt.query[,c("batch", "celltype", "celltype_proportion")],
  dt.atlas[,c("batch", "celltype", "celltype_proportion")]
)

##############################
## Dimensionality reduction ##
##############################

matrix <- dt %>% 
  dcast(batch~celltype, fill=0, value.var="celltype_proportion") %>% 
  matrix.please

# PCA
pca <- prcomp(matrix, rank.=5)

# UMAP
# umap <- uwot::umap(pca)

##########
## Plot ##
##########

to.plot <- pca$x %>% as.data.table %>% 
  .[,batch:=rownames(pca$x)] %>%
  merge(sample_metadata, by="batch")

# unique(to.plot$dataset)
stages <- unique(to.plot$stage) %>% sort
stage.colors <- viridis::viridis(n=length(stages))
names(stage.colors) <- rev(stages)
  
p <- ggplot(to.plot, aes(x=PC1, y=PC2, fill=stage, shape=dataset)) +
  geom_point(stroke=0.5, color="black", size=5) +
  scale_fill_manual(values=stage.colors) +
  scale_shape_manual(values=c(21,24)) +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  guides(shape=guide_legend(override.aes=list(fill="black"))) +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    axis.text = element_text(color="black", size=rel(0.8))
  )

pdf(paste0(io$outdir,"/pca_mapping_stages.pdf"), width=7, height=5, useDingbats = F)
print(p)
dev.off()

##########
## Save ##
##########

# fwrite(dt, io$outfile, sep="\t")

