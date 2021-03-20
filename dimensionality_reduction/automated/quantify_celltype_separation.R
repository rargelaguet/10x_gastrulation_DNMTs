ClusterPurity <- function(clusters, classes) {
  sum(apply(table(classes, clusters), 2, max)) / length(clusters)
}

suppressPackageStartupMessages(library(cluster))
# suppressPackageStartupMessages(library(clues))

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
} else {
  source("/homes/ricard/10x_gastrulation_DNMTs/settings.R")
}

# I/O
io$indir <- paste0(io$basedir,"/results/rna/dimensionality_reduction")
io$outdir <- paste0(io$basedir,"/results/rna/dimensionality_reduction/quantify_celltype_separation")

# Options
opts$samples <- list(
  # "E7.5_rep1" = "E7.5_rep1",
  # "E7.5_rep2" = "E7.5_rep2",
  # "E8.5_rep1" = "E8.5_rep1",
  # "E8.5_rep2" = "E8.5_rep2",
  # "E7.5" = c("E7.5_rep1","E7.5_rep2"),
  "E8.5" = c("E8.5_rep1","E8.5_rep2")
  # "all_cells" = c("E7.5_rep1","E7.5_rep2","E8.5_rep1","E8.5_rep2")
)

opts$npcs <- c(25,50)
opts$nfeatures <- c(1000,2000,3000)
opts$batch.correction <- "sample"


##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE] %>%
  .[,celltype.mapped:=as.factor(celltype.mapped)]

#################################
## Load latent representations ##
#################################

dt <- names(opts$samples) %>% map(function(i) {
  opts$nfeatures %>% map(function(j) {
    opts$npcs %>% map(function(k) {
      file <- sprintf("%s/%s/%s_pca_features%d_pcs%d.txt.gz",io$indir, i, paste(opts$samples[[i]],collapse="-"), j, k)
      # file <- sprintf("%s/%s_pca_features%d_pcs%d_batchcorrectionby%s.txt.gz",args$outdir, paste(args$samples,collapse="-"), args$features, args$npcs,paste(args$batch.correction,collapse="-"))
      if (file.exists(file)) {
        dt <- fread(file) %>% 
          merge(sample_metadata[,c("cell","celltype.mapped")], by="cell")  %>%
          .[,N:=.N,by="celltype.mapped"] %>% .[N>25] %>% .[,N:=NULL] %>% droplevels
        
        # calculate silhouette
        sh <- silhouette(as.numeric(dt$celltype.mapped), dist(dt[,grep("PC",colnames(dt)),with=F]))
        
        data.table(
          sample = i,
          nfeatures = j,
          npcs = k,
          silhouette = summary(sh)$avg.width %>% round(3)
        ) %>% return
      } else {
        print(sprintf("%s does not exist",file))
        return(NULL)
      }
    }) %>% rbindlist
  }) %>% rbindlist 
}) %>% rbindlist

# Save
# fwrite(dt, sprintf("%s_clustering.txt", args$outprefix), quote = TRUE)

##########
## Plot ##
##########

for (i in names(opts$samples)) {
  
  to.plot <- dt[sample==i] %>%
    .[,nfeatures:=factor(nfeatures,levels=opts$nfeatures)] %>%
    .[,npcs:=factor(npcs,levels=opts$npcs)] 
    
  # Heatmaps
  
  p <- ggplot(to.plot, aes(x=nfeatures, y=npcs)) +
    geom_tile(aes(fill=silhouette)) +
    labs(x="Number of features", y="Number of PCs") +
    scale_fill_gradient(low = "white", high = "red") +
    theme_classic() +
    theme(
      axis.text.y = element_text(colour="black",size=rel(1.0)),
      axis.text.x = element_text(colour="black",size=rel(1.0)),
      axis.ticks = element_blank()
    )
  
  pdf(sprintf("%s/%s_heatmap_silhouette.pdf",io$outdir,i))
  print(p)
  dev.off()
  
  # Bar plots
  
  to.plot2 <- to.plot %>%
    .[,group:=sprintf("%s_%s",nfeatures,npcs)]
  
  p <- ggplot(to.plot2, aes(x = nfeatures, y = silhouette, fill=npcs)) + 
    geom_bar(stat = "identity", position="dodge", color="black") +
    labs(x="Number of features", y="Silhouette") +
    theme_classic() +
    theme(
      axis.text.y = element_text(colour="black",size=rel(1.0)),
      axis.text.x = element_text(colour="black",size=rel(1.0)),
      axis.ticks.x = element_blank()
    )
    
  pdf(sprintf("%s/%s_barplots_silhouette.pdf",io$outdir,i))
  print(p)
  dev.off()
}
  
#################
## Save output ##
#################
