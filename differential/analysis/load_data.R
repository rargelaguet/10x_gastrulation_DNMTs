
##############
## Settings ##
##############

opts$celltypes = c(
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
  "Blood_progenitors",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid",
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
  "ExE_endoderm"
  # "ExE_ectoderm"
  # "Parietal_endoderm"
)

# Options
opts$ko.classes <- c(
  "Dnmt3a_KO", 
  # "Dnmt3a_HET_Dnmt3b_KO", 
  # "Dnmt3a_HET_Dnmt3b_WT", 
  # "Dnmt3a_KO_Dnmt3b_HET", 
  # "Dnmt3ab_KO", 
  "Dnmt3b_KO",
  "Dnmt1_KO"
)
opts$wt.class <- "WT"

# Define statistical significance
opts$threshold_fdr <- 0.01
opts$min.logFC <- 1

#############################################
## Load results from differential analysis ##
#############################################

# i <- opts$groupA[1]; j <- opts$celltypes[1]
diff.dt <- opts$ko.classes %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- sprintf("%s/%s/%s_%s_vs_%s.txt.gz", io$indir,i,j,opts$wt.class,i)
  if (file.exists(file)) {
    fread(file, select=c(1,2,4,6,7)) %>% .[,c("celltype","class"):=list(j,i)]
  }
}) %>% rbindlist }) %>% rbindlist %>%
  .[,celltype:=factor(celltype,levels=opts$celltypes)] %>%
  .[,class:=factor(class,levels=opts$ko.classes)] %>% 
  .[, sign := ifelse(logFC>0,sprintf("Upregulated in %s",class),sprintf("Downregulated in %s",class))] %>%
  .[, sig := (padj_fdr<=opts$threshold_fdr & abs(logFC)>=opts$min.logFC)]

# Print stats
print(sprintf("Number of classes: %s",length(unique(diff.dt$class))))
print(sprintf("Number of celltypes: %s",length(unique(diff.dt$celltype))))
print(sprintf("Number of genes: %s",length(unique(diff.dt$gene))))
