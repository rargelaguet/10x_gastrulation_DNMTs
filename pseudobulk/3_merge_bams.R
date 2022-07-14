  source(here::here("settings.R"))
  
  
  ######################
  ## Define arguments ##
  ######################
  
  p <- ArgumentParser(description='')
  
  
  # p$add_argument('--indir',      type="character",    help='Parent folder containing folders with barcode lists and param.rds files')
  # 
  # 
  # args <- p$parse_args(commandArgs(TRUE))
  args=list()
  args$indir <- ""
  
  bams <- dir(args$indir, pattern = ".bam$", full = TRUE)
  names <- basename(bams)
  outfiles <- unique(names)
  
  walk(outfiles, ~{
    outfile <- paste0(args$outdir, "/", .x)
    i <- names %in% .x
    infiles <- bams[i]
    
    cmd <- paste(
      "samtools merge",
      "-o", outfile,
      paste(infiles, collapse = " ")
    )
    print(cmd)
    system(cmd)
  })