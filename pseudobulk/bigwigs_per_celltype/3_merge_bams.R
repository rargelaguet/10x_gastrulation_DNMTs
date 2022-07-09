source(here::here("settings.R"))


# script to merge bam files of the same type
# alias file is used to match up files to be merged

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')

p$add_argument('--aliasfile',    type="character",    help='alias file')
args <- p$parse_args(commandArgs(TRUE))


###############
## Load data ##
###############

# Load alias file

tomerge <- fread(args$aliasfile) %>% 
  .[, .(sample, split_to, merge_to, outdir)] %>% 
  unique() %>% 
  split(by = "merge_to")




walk(tomerge, ~{
  outfile <- .x[, unique(merge_to)]
  soutfile <- gsub(".bam$", "_sorted.bam", outfile)
  
  if (file.exists(soutfile)) {
    print("already processed. skipping..")
    return(NULL)
  }
  
  
  infiles <- .x[, split_to]
  
  cmd1 <- paste(
    "samtools merge -f",
    outfile,
    paste(infiles, collapse = " ")
  )
  print(cmd1)
  system(cmd1)
  
  cmd2 <- paste(
    "samtools sort",
    outfile,
    ">", soutfile
  )
  print(cmd2)
  system(cmd2)
  
  # cmd3 <- paste(
  #   "samtools index",
  #   outfile
  # )
  # 
  # print(cmd3)
  # system(cmd3)
})
