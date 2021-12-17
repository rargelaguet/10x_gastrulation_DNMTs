#########
## I/O ##
#########

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
  io$script <- "/Users/ricard/10x_gastrulation_DNMTs/differential/differential.R"
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/10x_gastrulation_DNMTs/settings.R")
  io$script <- "/homes/ricard/10x_gastrulation_DNMTs/differential/differential.R"
  io$tmpdir <- "/hps/nobackup2/research/stegle/users/ricard/10x_gastrulation_DNMTs/results/differential/tmp"; dir.create(io$tmpdir, showWarnings=F)
} else {
  stop("Computer not recognised")
}
io$outdir <- paste0(io$basedir,"/results/differential"); dir.create(io$outdir, showWarnings = F)

#############
## Options ##
#############

# Testing mode
opts$test_mode <- FALSE

# Define groups
opts$query.groups <- c(
  # "E12.5_Dnmt3aWT_Dnmt3bHET",
  # "E12.5_Dnmt3aWT_Dnmt3bKO",
  # "E12.5_Dnmt3aHET_Dnmt3bWT",
  # "E12.5_Dnmt3aKO_Dnmt3bWT",
  "E8.5_Dnmt3aKO_Dnmt3bWT", 
  "E8.5_WT", 
  "E8.5_Dnmt3aHET_Dnmt3bKO", 
  "E8.5_Dnmt3aHET_Dnmt3bWT", 
  "E8.5_Dnmt3aKO_Dnmt3bHET", 
  "E8.5_Dnmt3aKO_Dnmt3bKO", 
  "E8.5_Dnmt3aWT_Dnmt3bKO",
  "E8.5_Dnmt1KO"
)

# Define query and reference group
opts$reference.groups <- "E8.5_WT"
opts$query.groups <- opts$classes[opts$classes!=opts$reference.groups]# %>% head(n=1)

# Define cell types with sufficient number of cells for differential comparison
opts$min.cells <- 50

opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  "PGC",
  "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  "Mixed_mesoderm",
  "Intermediate_mesoderm",
  "Caudal_Mesoderm",
  "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "Pharyngeal_mesoderm",
  "Cardiomyocytes",
  "Allantois",
  "ExE_mesoderm",
  "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm",
  "Parietal_endoderm"
)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  .[pass_QC==TRUE & class%in%opts$classes & celltype.mapped%in%opts$celltypes]

# Filter cell types with small number of cells in the reference group
opts$celltypes <- sample_metadata %>%
  .[class==opts$reference.groups,.N,by="celltype.mapped"] %>%
  .[N>=opts$min.cells,celltype.mapped]

#########
## Run ##
#########

for (i in opts$query.groups) {
  
  # Extract cell types
  opts$celltypes <- sample_metadata %>% 
    .[!is.na(celltype.mapped) & class%in%c(i,opts$reference.groups)] %>%
    .[,.N,by=c("class","celltype.mapped")] %>%
    .[N>=opts$min.cells] %>% .[,.N,by="celltype.mapped"] %>% .[N==2,celltype.mapped]
  
    for (j in opts$celltypes) {
      outfile <- sprintf("%s/%s_vs_%s_%s.txt.gz", io$outdir,i,opts$reference.groups,j)

      # Define LSF command
      if (grepl("ricard",Sys.info()['nodename'])) {
        lsf <- ""
      } else if (grepl("ebi",Sys.info()['nodename'])) {
        lsf <- sprintf("bsub -M 15000 -n 1 -q research-rh74 -o %s/%s_vs_%s_%s.txt", io$tmpdir,i,opts$reference.groups,j)
      }
      cmd <- sprintf("%s Rscript %s --groupA %s --groupB %s --celltype %s --outfile %s", lsf, io$script, opts$reference.groups, i, j, outfile)
      if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test_mode")

      # Run
      print(cmd)
      system(cmd)
    }
}

