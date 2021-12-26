here::i_am("processing/5_parse_sample_metadata_after_doublets.R")
source(here::here("settings.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',         type="character",   help='Metadata file to use as input')
p$add_argument('--doublet_files',    type="character", nargs="+",  help='Results of the doublet score detection algorithm')
p$add_argument('--outdir',          type="character",   help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

###################
## Load settings ##
###################

## START TEST ##
# args$metadata <- file.path(io$basedir,"results/qc/sample_metadata_after_qc.txt.gz")
# args$doublet_files <- list.files(file.path(io$basedir,"results_new/doublet_detection"), pattern = "*.txt.gz", full.names = T)
# args$outdir <- file.path(io$basedir,"results_new/doublet_detection/sample_metadata_after_doublets.txt.gz")
## END TEST ##

##########################
## Load doublet results ##
##########################

doublet.dt <- args$doublet_files %>% 
  map(~ fread(., select=c("cell","doublet_score","doublet_call"))) %>% 
  rbindlist

########################
## Plot doublet stats ##
########################

# to.plot <- doublet.dt[,mean(doublet_call),by="sample"]
# 
# p <- ggbarplot(to.plot, x="sample", y="V1", fill="gray70") +
#   # scale_fill_manual(values=opts$stage.colors) +
#   labs(x="", y="Fraction of doublets") +
#   coord_cartesian(ylim=c(0,0.10)) +
#   # facet_wrap(~stage)
#   theme(
#     legend.position = "none",
#     axis.text.y = element_text(colour="black",size=rel(0.8)),
#     axis.text.x = element_text(colour="black",size=rel(0.50), angle=20, hjust=1, vjust=1),
#   )
# 
# pdf(file.path(args$outdir,"doublets_barplot.pdf"), width=8, height=6)
# print(p)
# dev.off()


####################
## Merge and save ##
####################

# to.save <- fread(args$metadata) %>% 
#   merge(doublet.dt[,c("cell","hybrid_score","doublet_call")] %>% setnames("hybrid_score","doublet_score"), by="cell", all.x=TRUE)
to.save <- fread(args$metadata) %>% merge(doublet.dt, by="cell", all.x=TRUE)
fwrite(to.save, file.path(args$outdir,"sample_metadata_after_doublets.txt.gz"), sep="\t", na="NA", quote=F)
