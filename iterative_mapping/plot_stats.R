
#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/10x_gastrulation_DNMTs/settings.R")
}  
io$mapping.dir <- paste0(io$basedir,"/results/iterative_mapping")
io$outdir <- paste0(io$basedir,"/results/iterative_mapping/pdf")

opts$batches <- c(
  # "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001",
  # "SIGAB6_E85_3_Dnmt3aWT_Dnmt3b_WT_L002",
  # "SIGAC6_E85_5_Dnmt3aKO_Dnmt3b_Het_L003",
  # "SIGAD6_E85_8_Dnmt3aHet_Dnmt3b_KO_L004",
  "15_E8_5_D3A_WT_D3B_WT_L007",
  # "17_E8_5_D3A_KO_D3B_WT_L008",
  # "2_E8_5_D3A_WT_D3B_KO_L003",
  # "3_E8_5_D3A_HET_D3B_WT_L004",
  # "7_E8_5_D3A_WT_D3B_KO_L005",
  "8_E8_5_D3A_KO_D3B_KO_L006"
)

# Update sample metadata
sample_metadata <- sample_metadata %>%
  .[batch%in%opts$batches]
table(sample_metadata$batch)

################
## Load data  ##
################

mapping_dt <- opts$batches %>% map(function(i) {
  rbind(
    fread(sprintf("%s/%s_standard_mnn.txt.gz",io$mapping.dir,i)) %>% .[,batch:=i] %>% .[,method:="Standard MNN"]
    # fread(sprintf("%s/%s_iterative_mnn.txt.gz",io$mapping.dir,i)) %>% .[,batch:=i] %>% .[,method:="Tree-guided MNN"]
  )
}) %>% rbindlist

mapping_dt %>%
  .[,celltype_mapped:=stringr::str_replace_all(celltype_mapped," ", "_")] %>%
  .[,celltype_mapped:=factor(celltype_mapped,levels=opts$celltypes)]

unique(mapping_dt$celltype_mapped)

##########
## Plot ##
##########

barplot.pub <- function(df, x, colors=NULL, xlim.max=NULL) {
  p <- ggplot(df, aes_string(x=x, y="N")) +
    scale_x_discrete(drop=FALSE) + 
    labs(y="Number of cells") +
    # theme_classic() +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=rel(1.3)),
      axis.title.x = element_text(color="black", size=rel(1.1)),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1.3), color="black"),
      axis.text.x = element_text(size=rel(1.1), color="black")
    )
  
  if (is.null(colors)) {
    p <- p + geom_bar(stat="identity", color="black")
  } else {
    p <- p + geom_bar(aes_string(fill=x), stat="identity", color="black") + 
      scale_fill_manual(values=colors, drop=F)
  }
  
  if (!is.null(xlim.max)) {
    p <- p + coord_flip(ylim=c(0,xlim.max))
  } else {
    p <- p + coord_flip()
  }
  
  return(p)
}

for (i in opts$batches) {
  
  to.plot <- mapping_dt[batch==i] %>% 
    merge(sample_metadata, by=c("cell","batch")) %>% 
    .[celltype_mapped=="Forebrain Midbrain Hindbrain",celltype_mapped:="Forebrain_Midbrain_Hindbrain"] %>%
    .[,.N,by=c("celltype_mapped","batch","method")]
  
  p <- barplot.pub(to.plot, x="celltype_mapped", colors=opts$celltype.colors) +
    facet_wrap(~method, nrow=1, scales="free_x") +
    theme(
      strip.background = element_blank()
    )  
  
  # pdf(sprintf("%s/barplots_%s.pdf",io$outdir,i), width=18, height=7)
  pdf(sprintf("%s/barplots_%s.pdf",io$outdir,i), width=9, height=7)
  print(p)
  dev.off()
}
