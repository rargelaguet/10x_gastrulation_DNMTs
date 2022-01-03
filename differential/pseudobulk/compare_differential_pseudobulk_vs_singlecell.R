#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
} else {
  stop("Computer not recognised")
}


# I/O
io$outdir <- paste0(io$basedir,"/results/rna/differential/pseudobulk/pdf")
io$rna_diff_singlecell <- paste0(io$basedir,"/results/rna/differential")
io$rna_diff_pseudobulk <- paste0(io$basedir,"/results/rna/differential/pseudobulk")

# Options
opts$celltypes <- c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  # "PGC",
  "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  "Mixed_mesoderm",
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
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm",
  "Parietal_endoderm"
) %>% head(n=5)

#######################
## Load marker genes ##
#######################

marker_genes.dt <- fread(io$rna.atlas.marker_genes.up)
marker_genes <- unique(marker_genes.dt$gene)

###########################################################
## Load differential RNA expression at single-cell level ##
###########################################################

rna_diff_singlecell.dt <- opts$celltypes %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- sprintf("%s/%s_vs_%s.txt.gz", io$rna_diff_singlecell,i,j)
  if (file.exists(file)) {
    fread(file, select = c(1,2)) %>% 
      .[gene%in%marker_genes] %>%
      # .[is.na(logFC),logFC:=0] %>%
      .[,c("celltypeA","celltypeB"):=list(as.factor(i),as.factor(j))] %>%
      return
  }
}) %>% rbindlist }) %>% rbindlist

# %>%
#   .[,celltypeA:=factor(celltypeA,levels=opts$celltypes)] %>%
#   .[,celltypeB:=factor(celltypeB,levels=opts$celltypes)]# %>%
  # .[,sig:=abs(logFC)>=opts$min.Log2FC & padj_fdr<=opts$min.FDR] %>%
  # .[is.na(sig),sig:=FALSE]


##########################################################
## Load differential RNA expression at pseudobulk level ##
##########################################################

rna_diff_pseudobulk.dt <- opts$celltypes %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- sprintf("%s/%s_vs_%s.txt.gz", io$rna_diff_pseudobulk,i,j)
  if (file.exists(file)) {
    fread(file, select = c(1,2,3)) %>% 
      .[gene%in%marker_genes] %>%
      .[,logFC:=round(groupB-groupA,2)] %>% .[,c("groupA","groupB"):=NULL] %>%
      .[,c("celltypeA","celltypeB"):=list(as.factor(i),as.factor(j))] %>%
      return
  }
}) %>% rbindlist }) %>% rbindlist

#############
## Combine ##
#############

rna_diff.dt <- rbind(
  rna_diff_singlecell.dt[,class:=as.factor("singlecell")],
  rna_diff_pseudobulk.dt[,class:=as.factor("pseudobulk")]
) %>% dcast(gene+celltypeA+celltypeB~class, value.var="logFC") %>%
  .[!is.na(pseudobulk) & !is.na(singlecell)]

##########
## Plot ##
##########

opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors)%in%opts$celltypes]


for (i in 1:length(opts$celltypes)) {
  for (j in i:length(opts$celltypes)) {
    if (i!=j) {
      
      to.plot <- rna_diff.dt[celltypeA==opts$celltypes[i] & celltypeB==opts$celltypes[j]]
            
      p <- ggscatter(to.plot, x="singlecell", y="pseudobulk", fill="gray80", shape=21, size=2, 
                      add="reg.line", add.params = list(color="black", fill="lightgray", alpha=0.5, size=.5), conf.int=TRUE) +
        stat_cor(method = "pearson", label.x.npc = "middle", label.y.npc = "bottom") +
        labs(x=sprintf("%s RNA DE (single-cell)",opts$celltypes[i]), y=sprintf("%s RNA DE (pseudobulk)",opts$celltypes[j])) +
        theme(
          legend.position = "none",
          axis.text = element_text(size=rel(0.7)),
          axis.title = element_text(size=rel(0.85))
        )
      
      # pdf(sprintf("%s/%s_chromvar_singlecell_vs_pseudobulk.pdf",io$outdir,i), width=13, height=5)
      print(p)
      # dev.off()
    }
  }
}
