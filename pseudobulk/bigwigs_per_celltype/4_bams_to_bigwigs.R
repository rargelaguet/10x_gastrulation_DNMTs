source(here::here("settings.R"))

library(Rsamtools)
library(rtracklayer)
library(GenomicRanges)
library(GenomicAlignments)
library(BSgenome.Mmusculus.UCSC.mm10)

######################
## Define arguments ##
######################

#p <- ArgumentParser(description='')


p$add_argument('--indir',      type="character",    help='Folder containing bam files')
p$add_argument('--outdir',      type="character",    help='Folder to save bigwig files')


args <- p$parse_args(commandArgs(TRUE))


bams <- dir(args$indir, pattern = "_sorted.bam$", full = TRUE)


bigwigs <- paste0(
  args$outdir,
  "/",
  gsub("_sorted.bam$", ".bw", basename(bams))
)


walk2(bams, bigwigs, ~{
  print(basename(.x))
  if (file.exists(.y)){
    print("already processed. skipping..")
    return(NULL)
  }
  
  readBam <- possibly(readGAlignments, otherwise = NA)
  gr <- readBam(.x)
  if (length(gr)<2) {
    print("Bam file error. Skipping..")
    return(NULL)
  }
  
  gr <- gr %>% 
    as("GRanges") %>% 
    resize(90) %>%  # resize so all reads are the 90bp 10x read length
    trim()
  
  libsize <- length(gr)
  
  cov <- coverage(gr)
  
  rpm <- map(as.list(cov), ~signif(10^6 * .x/libsize, digits = 3)) %>% 
    as("SimpleRleList")
  
  
  export.bw(rpm, .y)
})