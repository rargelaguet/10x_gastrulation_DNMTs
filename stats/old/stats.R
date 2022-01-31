here::i_am("processing/stats/stats.R")

source(here::here("settings.R"))

#####################
## Define settings ##
#####################

# I/O
# Options
opts$classes <- c(
    "WT",
    "Dnmt3a_KO",
    "Dnmt3b_KO",
    "Dnmt1_KO",
    "Dnmt3ab_KO"
)

###############
## Load data ##
###############

metadata <- fread(io$metadata) %>% 
    .[,dataset:=ifelse(grepl("Grosswendt",sample),"CRISPR (Grosswendt2020)","KO")] %>%
    .[pass_rnaQC==TRUE & class%in%opts$classes]

#################################################
## Plot number of cells per class and data set ##
#################################################

to.plot <- metadata %>%
    .[,.N,by=c("class","dataset")]

p <- ggbarplot(to.plot, x="class", y="N", fill="dataset", position=position_dodge(width = 0.75)) +
    # scale_fill_manual(values=opts$stage.colors) +
    labs(x="", y="Number of cells") +
    theme(
        legend.position = "right",
        axis.text.y = element_text(colour="black",size=rel(0.8)),
        axis.text.x = element_text(colour="black",size=rel(0.80)),
    )

pdf(file.path(args$outdir,"qc_metrics_barplot.pdf"), width=8, height=6)
print(p)
dev.off()

##########
## Save ##
##########

fwrite(metadata, file.path(args$outdir,"sample_metadata_after_qc.txt.gz"), quote=F, na="NA", sep="\t")

