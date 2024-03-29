here::i_am("differential/run.R")

source(here::here("settings.R"))

#####################
## Define settings ##
#####################

io$script <- here::here("differential/differential.R")
io$outdir <- file.path(io$basedir,"results_new/differential/dnmt1_vs_dnmt3ab_ko"); dir.create(io$outdir, showWarnings=F)

# Rename celltypes
opts$rename_celltypes <- c(
  "Erythroid3" = "Erythroid",
  "Erythroid2" = "Erythroid",
  "Erythroid1" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors"
  # "Intermediate_mesoderm" = "Mixed_mesoderm",
  # "Paraxial_mesoderm" = "Mixed_mesoderm",
  # "Nascent_mesoderm" = "Mixed_mesoderm",
  # "Pharyngeal_mesoderm" = "Mixed_mesoderm"
  # "Visceral_endoderm" = "ExE_endoderm"
)

opts$ko.classes <- c(
  "E8.5_Dnmt3aKO_Dnmt3bKO"
)

opts$wt.class <- "E8.5_Dnmt1KO"

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,opts$rename_celltypes)] %>%
  .[pass_rnaQC==TRUE & !is.na(celltype.mapped)]

# Only consider cell types with sufficient observations in WT cells
celltypes.to.use <- sample_metadata[class==opts$wt.class,.(N=.N),by="celltype.mapped"] %>% .[N>=50,celltype.mapped]
sample_metadata <- sample_metadata[celltype.mapped%in%celltypes.to.use]

print(table(sample_metadata$celltype.mapped,sample_metadata$class))

###################################
## Run all pair-wise comparisons ##
###################################

opts$min.cells <- 30

# j <- "Blood_progenitors"; i <- "E8.5_Dnmt3aKO_Dnmt3bKO"
for (i in opts$ko.classes) {
  
  # Define cell types to use 
  celltypes.to.use <- sample_metadata %>% .[class==i,.N,by="celltype.mapped"] %>% .[N>=opts$min.cells,celltype.mapped]
  
  for (j in celltypes.to.use) {
    outfile <- sprintf("%s/%s_%s_vs_%s.txt.gz", io$outdir,j,opts$wt.class,i); dir.create(dirname(outfile), showWarnings = F)
    if (!file.exists(outfile)) {
      
      # Define LSF command
      if (grepl("BI",Sys.info()['nodename'])) {
        lsf <- ""
      } else if (grepl("pebble|headstone", Sys.info()['nodename'])) {
        lsf <- sprintf("sbatch -n 1 --mem 30G --wrap")
      }
      cmd <- sprintf("%s 'Rscript %s --groupA %s --groupB %s --celltype %s --group_label class --outfile %s'", lsf, io$script, opts$wt.class, i, j, outfile)
      
      # Run
      print(cmd)
      system(cmd)
    }
  }
}

