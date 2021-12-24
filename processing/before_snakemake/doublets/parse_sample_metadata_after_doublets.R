###################
## Load settings ##
###################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
} else {
  source("/homes/ricard/10x_gastrulation_DNMTs/settings.R")
}

#########
## I/O ##
#########

io$metadata <- paste0(io$basedir,"/results/mapping/sample_metadata_after_mapping.txt.gz")
io$doublet.dir <- paste0(io$basedir,"/results/doublets")
io$output.metadata <- paste0(io$doublet.dir,"/sample_metadata_after_doublets.txt.gz")

#############
## Options ##
#############

# opts$samples <- c(
# 	"E7.5_rep1",
# 	"E7.5_rep2", 
# 	"E8.5_rep1",
# 	"E8.5_rep2"
# )

opts$hybrid_score_threshold <- 1

###############
## Load data ##
###############

# Load doublet calls
doublets.dt <- opts$batches %>% 
  map(~ fread(sprintf("%s/%s_%s.txt.gz",io$doublet.dir,.,opts$hybrid_score_threshold))) %>% 
  rbindlist

# Merge with sample metadata ##
sample_metadata <- fread(io$metadata) %>% 
  merge(doublets.dt,by=c("cell","batch"), all.x=TRUE)

head(sample_metadata)

#################
## Save output ##
#################

fwrite(sample_metadata, io$output.metadata, sep="\t", na="NA", quote=F)
