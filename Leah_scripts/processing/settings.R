io <- list()

io$samples <- list()

io$samples$second_batch <- c(
  "E125_DNMT3A_HET_A_L001",
  "E125_DNMT3A_HET_A_L003",
  "E125_DNMT3A_KO_B_L002",
  "E125_DNMT3A_KO_E_L004"
)

io$samples$third_batch <- c(
  "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001",
  "SIGAB6_E85_3_Dnmt3aWT_Dnmt3b_WT_L002",
  "SIGAC6_E85_5_Dnmt3aKO_Dnmt3b_Het_L003",
  "SIGAD6_E85_8_Dnmt3aHet_Dnmt3b_KO_L004"
)

io$samples$fourth_batch <- c(
  "2_E8_5_D3A_WT_D3B_KO_L003",
  "3_E8_5_D3A_HET_D3B_WT_L004",
  "7_E8_5_D3A_WT_D3B_KO_L005",
  "8_E8_5_D3A_KO_D3B_KO_L006",
  "15_E8_5_D3A_WT_D3B_WT_L007",
  "17_E8_5_D3A_KO_D3B_WT_L008",
  "A_E12_5_D3a_Het_L001",
  "B_E12_5_D3a_KO_L002"
)

io$subset.proteincoding <- NULL # Choose between NULL and path to proteincoding genes
io$qc <- TRUE # Choose between NULL and TRUE

#################
## QC Settings ##
#################

io$min_nFeature_RNA <- list()
io$min_nCount_RNA <- list()
io$max_percent.mt <- list()
io$max_doublet_score <- list()

##################
## Second Batch ##
##################

io$min_nFeature_RNA$second_batch <- c(
  "E125_DNMT3A_HET_A_L001" = 0,
  "E125_DNMT3A_HET_A_L003" = 0,
  "E125_DNMT3A_KO_B_L002" = 1500,
  "E125_DNMT3A_KO_E_L004" = 500
)
io$min_nCount_RNA$second_batch <- c(
  "E125_DNMT3A_HET_A_L001" = 200,
  "E125_DNMT3A_HET_A_L003" = 200,
  "E125_DNMT3A_KO_B_L002" = 2000,
  "E125_DNMT3A_KO_E_L004" = 1000
)
io$max_percent.mt$second_batch <- c(
  "E125_DNMT3A_HET_A_L001" = 50,
  "E125_DNMT3A_HET_A_L003" = 50,
  "E125_DNMT3A_KO_B_L002" = 25,
  "E125_DNMT3A_KO_E_L004" = 25
)
io$max_doublet_score$second_batch <- c(
  "E125_DNMT3A_HET_A_L001" = 10000,
  "E125_DNMT3A_HET_A_L003" = 3000,
  "E125_DNMT3A_KO_B_L002" = 2000,
  "E125_DNMT3A_KO_E_L004" = 10000
)

#################
## Third Batch ##
#################

io$min_nFeature_RNA$third_batch <- c(
  "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001" = 3000,
  "SIGAB6_E85_3_Dnmt3aWT_Dnmt3b_WT_L002" = 3000,
  "SIGAC6_E85_5_Dnmt3aKO_Dnmt3b_Het_L003" = 3500,
  "SIGAD6_E85_8_Dnmt3aHet_Dnmt3b_KO_L004" = 3500
)
io$min_nCount_RNA$third_batch <- c(
  "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001" = 6000,
  "SIGAB6_E85_3_Dnmt3aWT_Dnmt3b_WT_L002" = 5000,
  "SIGAC6_E85_5_Dnmt3aKO_Dnmt3b_Het_L003" = 8500,
  "SIGAD6_E85_8_Dnmt3aHet_Dnmt3b_KO_L004" = 5000
)
io$max_percent.mt$third_batch <- c(
  "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001" = 10,
  "SIGAB6_E85_3_Dnmt3aWT_Dnmt3b_WT_L002" = 10,
  "SIGAC6_E85_5_Dnmt3aKO_Dnmt3b_Het_L003" = 10,
  "SIGAD6_E85_8_Dnmt3aHet_Dnmt3b_KO_L004" = 10
)
io$max_doublet_score$third_batch <- c(
  "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001" = 1000,
  "SIGAB6_E85_3_Dnmt3aWT_Dnmt3b_WT_L002" = 1000,
  "SIGAC6_E85_5_Dnmt3aKO_Dnmt3b_Het_L003" = 1000,
  "SIGAD6_E85_8_Dnmt3aHet_Dnmt3b_KO_L004" = 1000
)

##################
## Fourth Batch ##
##################

#io$min_nFeature_RNA$fourth_batch <- c(
#  "2_E8_5_D3A_WT_D3B_KO_L003" = 1215,
#  "3_E8_5_D3A_HET_D3B_WT_L004" = 1215,
#  "7_E8_5_D3A_WT_D3B_KO_L005" = 1215,
#  "8_E8_5_D3A_KO_D3B_KO_L006" = 1215,
#  "15_E8_5_D3A_WT_D3B_WT_L007" = 1215,
#  "17_E8_5_D3A_KO_D3B_WT_L008" = 1215,
#  "A_E12_5_D3a_Het_L001" = 1215,
#  "B_E12_5_D3a_KO_L002" = 1215
#)
io$min_nFeature_RNA$fourth_batch <- c(
  "2_E8_5_D3A_WT_D3B_KO_L003" = 2500,
  "3_E8_5_D3A_HET_D3B_WT_L004" = 3500,
  "7_E8_5_D3A_WT_D3B_KO_L005" = 2500,
  "8_E8_5_D3A_KO_D3B_KO_L006" = 3500,
  "15_E8_5_D3A_WT_D3B_WT_L007" = 3500,
  "17_E8_5_D3A_KO_D3B_WT_L008" = 3500,
  "A_E12_5_D3a_Het_L001" = 2500,
  "B_E12_5_D3a_KO_L002" = 3000
)
#io$min_nCount_RNA$fourth_batch <- c(
#  "2_E8_5_D3A_WT_D3B_KO_L003" = 5000,
#  "3_E8_5_D3A_HET_D3B_WT_L004" = 5000,
#  "7_E8_5_D3A_WT_D3B_KO_L005" = 5000,
#  "8_E8_5_D3A_KO_D3B_KO_L006" = 5000,
#  "15_E8_5_D3A_WT_D3B_WT_L007" = 5000,
#  "17_E8_5_D3A_KO_D3B_WT_L008" = 5000,
#  "A_E12_5_D3a_Het_L001" = 5000,
#  "B_E12_5_D3a_KO_L002" = 5000
#)
io$min_nCount_RNA$fourth_batch <- c(
  "2_E8_5_D3A_WT_D3B_KO_L003" = 6000,
  "3_E8_5_D3A_HET_D3B_WT_L004" = 10000,
  "7_E8_5_D3A_WT_D3B_KO_L005" = 5000,
  "8_E8_5_D3A_KO_D3B_KO_L006" = 10000,
  "15_E8_5_D3A_WT_D3B_WT_L007" = 10000,
  "17_E8_5_D3A_KO_D3B_WT_L008" = 10000,
  "A_E12_5_D3a_Het_L001" = 10000,
  "B_E12_5_D3a_KO_L002" = 6000
)
io$max_percent.mt$fourth_batch <- c(
  "2_E8_5_D3A_WT_D3B_KO_L003" = 10,
  "3_E8_5_D3A_HET_D3B_WT_L004" = 10,
  "7_E8_5_D3A_WT_D3B_KO_L005" = 10,
  "8_E8_5_D3A_KO_D3B_KO_L006" = 10,
  "15_E8_5_D3A_WT_D3B_WT_L007" = 10,
  "17_E8_5_D3A_KO_D3B_WT_L008" = 10,
  "A_E12_5_D3a_Het_L001" = 10,
  "B_E12_5_D3a_KO_L002" = 10
)
io$max_doublet_score$fourth_batch <- c(
  "2_E8_5_D3A_WT_D3B_KO_L003" = 1000,
  "3_E8_5_D3A_HET_D3B_WT_L004" = 1000,
  "7_E8_5_D3A_WT_D3B_KO_L005" = 1000,
  "8_E8_5_D3A_KO_D3B_KO_L006" = 1000,
  "15_E8_5_D3A_WT_D3B_WT_L007" = 1000,
  "17_E8_5_D3A_KO_D3B_WT_L008" = 1000,
  "A_E12_5_D3a_Het_L001" = 1000,
  "B_E12_5_D3a_KO_L002" = 1000
)