matrix.please <- function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[[1]]
  m
}

#####################
## Define settings ##
#####################

source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
io$outdir <- paste0(io$basedir,"/results/mapping_stages")

sample_metadata <- sample_metadata[stage!="E12.5"]

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
dt.atlas <- fread(paste0(io$atlas.basedir,"/results/general_stats/celltype_proportions.txt.gz")) %>%
  .[sample%in%sample_metadata.atlas$batch] %>%
  setnames("sample","batch")

#################
## Concatenate ##
#################

sample_metadata <- rbind(sample_metadata.query, sample_metadata.atlas)
dt <- rbind(dt.query,dt.atlas)

##############################
## Dimensionality reduction ##
##############################

matrix <- dt %>% dcast(batch~celltype, fill=0, value.var="celltype_proportion") %>% 
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

