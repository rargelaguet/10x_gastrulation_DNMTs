#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
} else {
  source("/homes/ricard/10x_gastrulation_DNMTs/settings.R")
}


#############
## Options ##
#############

opts$batches <- c(
  # E12.5  
  # "E125_DNMT3A_HET_A_L001",
  # "E125_DNMT3A_HET_A_L003",
  # "E125_DNMT3A_KO_B_L002",
  # "E125_DNMT3A_KO_E_L004",
  # "A_E12_5_D3a_Het_L001",
  # "B_E12_5_D3a_KO_L002",
  
  
  # E8.5  
  # "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001",
  # "SIGAB6_E85_3_Dnmt3aWT_Dnmt3b_WT_L002",
  # "SIGAC6_E85_5_Dnmt3aKO_Dnmt3b_Het_L003",
  # "SIGAD6_E85_8_Dnmt3aHet_Dnmt3b_KO_L004",
  # "15_E8_5_D3A_WT_D3B_WT_L007",
  # "17_E8_5_D3A_KO_D3B_WT_L008",
  # "2_E8_5_D3A_WT_D3B_KO_L003",
  # "3_E8_5_D3A_HET_D3B_WT_L004",
  # "7_E8_5_D3A_WT_D3B_KO_L005",
  # "8_E8_5_D3A_KO_D3B_KO_L006",
  # "E8_5_Dnmt1_KO_male_SIGAC8_L001",
  # "E8_5_Dnmt1_KO_male_SIGAD8_L002",
  # "E8_5_Dnmt1_KO_male_SIGAE8_L003",
  # "E8_5_Dnmt1_WT_female_SIGAB8_L004",
  # "E8_5_Dnmt1_WT_female_SIGAF8_L005",
  # "E8_5_Dnmt3ab_WT_female_SIGAA8_L006",
  # "SIGAH10_Dnmt3ab_WT_L002",
  # "SIGAH11_Dnmt3ab_WT_L003",
  # "SIGAH9_Dnmt3a_KO_Dnmt3b_Het_L001"
  "SIGAG5_9_dnmt3ab_DKO_L005"
)

#########
## I/O ##
#########

io$metadata <- "/hps/nobackup2/research/stegle/users/ricard/10x_gastrulation_DNMTs/processed/sixth_batch/metadata.txt.gz"

io$mapping.dir <- paste0(io$basedir,"/results/mapping")
io$outfile <- io$metadata

###############
## Load data ##
###############

# Load cell metadata
sample_metadata <- fread(io$metadata)# %>% .[batch%in%opts$batches]

# Load mapping results
mapping.dt <- opts$batches %>% map(function(i) {
  readRDS(sprintf("%s/mapping_mnn_%s.rds",io$mapping.dir,i))$mapping %>% 
    .[,c("cell","celltype.mapped")] %>% as.data.table %>% 
    .[,batch:=i]
    # .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,"_"," ")]
}) %>% rbindlist

# bar <- fread("/Users/ricard/data/10x_gastrulation_DNMTs/results/second_batch/mapping/leah/mapping_mnn_WT.txt.gz") %>%
#   .[cell %in% mapping.dt$cell] %>% .[,c("cell","celltype.mapped")] %>%
#   .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,"/"," ")]
# foobar <- merge(foo,bar,by="cell")
# mean(foobar$celltype.mapped.x == foobar$celltype.mapped.y)
# disagreement <- foobar[!celltype.mapped.x == celltype.mapped.y]
# fwrite(disagreement, "/Users/ricard/data/10x_gastrulation_DNMTs/results/second_batch/mapping/disagreement_leah_vs_ricard.txt.gz")

###########
## Merge ##
###########

sample.metadata <- sample_metadata %>%
  merge(mapping.dt, by=c("cell","batch"), all.x=TRUE)

#################
## Save output ##
#################

fwrite(sample.metadata, io$outfile, sep="\t", na="NA", quote=F)


