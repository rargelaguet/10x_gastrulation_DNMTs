source(here::here("settings.R"))



# this script uses a series of samtools commands to split out a 10x 3' RNA bam file by barcode
# the barcode info is taken from the alias file 
# script needs to be run once for each 10x bam file
 
# samtools comands are taken from: https://kb.10xgenomics.com/hc/en-us/articles/360022448251-How-to-filter-the-BAM-file-produced-by-10x-pipelines-with-a-list-of-barcodes-

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')

p$add_argument('--aliasfile',    type="character",    help='alias file')
p$add_argument('--sample',    type="character",    help='10x sample name to be processed')


args <- p$parse_args(commandArgs(TRUE))



#####################
## Define settings ##
#####################



###############
## Load data ##
###############

# Load alias file

aliases <- fread(args$aliasfile)





sub_by_sample <- aliases[sample == args$sample]
types <- sub_by_sample[, unique(split_to)]

# iterate over celltypes
walk(types, ~{
  outfile <- .x
  print(paste("processing ", .x))
  sub_by_type <- sub_by_sample[celltype_class == .x]
  
  
  
  bam <- sub_by_type[, unique(bam_path)]
  
  outdir <- dirname(.x)
  tempdir <- paste0(outdir, "/temp")
  dir.create(tempdir, recursive = TRUE)
  
  
  
  if (file.exists(outfile)) {
    print("bam file already exists. skipping..")
    return(NULL)
  }
  
  
  
  barcodes <- sub_by_type[, .(barcode)]
  bcfile <- paste0(tempdir, "/barcodes.txt")
  fwrite(barcodes, bcfile, quote = FALSE, col.names = FALSE)
  
  headfile <- paste0(tempdir, "/header")
  filtfile <- paste0(tempdir, "/filt")
  samfile <- paste0(tempdir, "/filtered.sam")
  
  # Save the header lines
  cmd1 <- paste(
    "samtools view -H",
    bam,
    ">",
    headfile
  )
  print(cmd1)
  system(cmd1)
  
  # Filter alignments using filter.txt. Use LC_ALL=C to set C locale instead of UTF-8
  cmd2 <- paste(
    "samtools view",
    bam,
    "| LC_ALL=C grep -F -f",
    bcfile,
    ">",
    filtfile
    
  )
  print(cmd2)
  system(cmd2)
  
  
  # Combine header and body
  
  cmd3 <- paste(
    "cat",
    headfile,
    filtfile,
    ">",
    samfile
  )
  
  print(cmd3)
  system(cmd3)
  
  # Convert filtered.sam to BAM format
  
  cmd4 <- paste(
    "samtools view -b",
    samfile,
    ">",
    outfile
  )
  
  print(cmd4)
  system(cmd4)
  
  # cleanup
  cmd5 <- paste(
    "rm",
    bcfile,
    headfile,
    filtfile,
    samfile
  )
  
  print(cmd5)
  system(cmd5)
  
})

