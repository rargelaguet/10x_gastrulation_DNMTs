io <- list()

io$samples <- list()

io$samples$second_batch <- c(
  "E75_TET_TKO_L002",
  "E75_WT_Host_L001",
  "E85_Rep1_TET_TKO_L004",
  "E85_Rep2_TET_TKO_L006",
  "E85_Rep1_WT_Host_L003",
  "E85_Rep2_WT_Host_L005",
  "E125_DNMT3A_HET_A_L001",
  "E125_DNMT3A_HET_A_L003",
  "E125_DNMT3A_KO_B_L002",
  "E125_DNMT3A_KO_E_L004"
)

io$samples$third_batch <- c(
  "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001",
  "SIGAB6_E85_3_Dnmt3aWT_Dnmt3b_WT_L002",
  "SIGAC6_E85_5_Dnmt3aKO_Dnmt3b_Het_L003",
  "SIGAD6_E85_8_Dnmt3aHet_Dnmt3b_KO_L004",
  "SIGAE4_E105_3_TET123_Chimera_Host_L005",
  "SIGAF4_E105_3_TET123_Chimera_TKO_L006",
  "SIGAG4_E105_5_TET123_Chimera_Host_L007",
  "SIGAH4_E105_5_TET123_Chimera_TKO_L008"
)

io$subset.proteincoding <- NULL # Choose between NULL and path to proteincoding genes
io$qc <- TRUE # Choose between NULL and TRUE

io$min_nFeature_RNA <- list()
io$min_nCount_RNA <- list()
io$max_percent.mt <- list()
io$max_doublet_score <- list()
io$min_nFeature_RNA$second_batch <- c(
  "E75_TET_TKO_L002" = 3300,
  "E75_WT_Host_L001" = 2000,
  "E85_Rep1_TET_TKO_L004" = 2500,
  "E85_Rep2_TET_TKO_L006" = 2000,
  "E85_Rep1_WT_Host_L003" = 2000,
  "E85_Rep2_WT_Host_L005" = 2000,
  "E125_DNMT3A_HET_A_L001" = 0,
  "E125_DNMT3A_HET_A_L003" = 0,
  "E125_DNMT3A_KO_B_L002" = 1500,
  "E125_DNMT3A_KO_E_L004" = 500
)
io$min_nCount_RNA$second_batch <- c(
  "E75_TET_TKO_L002" = 14000,
  "E75_WT_Host_L001" = 5000,
  "E85_Rep1_TET_TKO_L004" = 6000,
  "E85_Rep2_TET_TKO_L006" = 5000,
  "E85_Rep1_WT_Host_L003" = 5000,
  "E85_Rep2_WT_Host_L005" = 5000,
  "E125_DNMT3A_HET_A_L001" = 200,
  "E125_DNMT3A_HET_A_L003" = 200,
  "E125_DNMT3A_KO_B_L002" = 2000,
  "E125_DNMT3A_KO_E_L004" = 1000
)
io$max_percent.mt$second_batch <- c(
  "E75_TET_TKO_L002" = 10,
  "E75_WT_Host_L001" = 10,
  "E85_Rep1_TET_TKO_L004" = 10,
  "E85_Rep2_TET_TKO_L006" = 10,
  "E85_Rep1_WT_Host_L003" = 10,
  "E85_Rep2_WT_Host_L005" = 10,
  "E125_DNMT3A_HET_A_L001" = 50,
  "E125_DNMT3A_HET_A_L003" = 50,
  "E125_DNMT3A_KO_B_L002" = 25,
  "E125_DNMT3A_KO_E_L004" = 25
)
io$max_doublet_score$second_batch <- c(
  "E75_TET_TKO_L002" = 1000,
  "E75_WT_Host_L001" = 1000,
  "E85_Rep1_TET_TKO_L004" = 1000,
  "E85_Rep2_TET_TKO_L006" = 1000,
  "E85_Rep1_WT_Host_L003" = 1000,
  "E85_Rep2_WT_Host_L005" = 1000,
  "E125_DNMT3A_HET_A_L001" = 10000,
  "E125_DNMT3A_HET_A_L003" = 3000,
  "E125_DNMT3A_KO_B_L002" = 2000,
  "E125_DNMT3A_KO_E_L004" = 10000
)

io$min_nFeature_RNA$third_batch <- c(
  "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001" = 3000,
  "SIGAB6_E85_3_Dnmt3aWT_Dnmt3b_WT_L002" = 3000,
  "SIGAC6_E85_5_Dnmt3aKO_Dnmt3b_Het_L003" = 3500,
  "SIGAD6_E85_8_Dnmt3aHet_Dnmt3b_KO_L004" = 3500,
  "SIGAE4_E105_3_TET123_Chimera_Host_L005" = 3500,
  "SIGAF4_E105_3_TET123_Chimera_TKO_L006" = 2000,
  "SIGAG4_E105_5_TET123_Chimera_Host_L007" = 3500,
  "SIGAH4_E105_5_TET123_Chimera_TKO_L008" = 2000
)
io$min_nCount_RNA$third_batch <- c(
  "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001" = 6000,
  "SIGAB6_E85_3_Dnmt3aWT_Dnmt3b_WT_L002" = 5000,
  "SIGAC6_E85_5_Dnmt3aKO_Dnmt3b_Het_L003" = 8500,
  "SIGAD6_E85_8_Dnmt3aHet_Dnmt3b_KO_L004" = 5000,
  "SIGAE4_E105_3_TET123_Chimera_Host_L005" = 10000,
  "SIGAF4_E105_3_TET123_Chimera_TKO_L006" = 8000,
  "SIGAG4_E105_5_TET123_Chimera_Host_L007" = 8000,
  "SIGAH4_E105_5_TET123_Chimera_TKO_L008" = 8000
)
io$max_percent.mt$third_batch <- c(
  "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001" = 10,
  "SIGAB6_E85_3_Dnmt3aWT_Dnmt3b_WT_L002" = 10,
  "SIGAC6_E85_5_Dnmt3aKO_Dnmt3b_Het_L003" = 10,
  "SIGAD6_E85_8_Dnmt3aHet_Dnmt3b_KO_L004" = 10,
  "SIGAE4_E105_3_TET123_Chimera_Host_L005" = 10,
  "SIGAF4_E105_3_TET123_Chimera_TKO_L006" = 10,
  "SIGAG4_E105_5_TET123_Chimera_Host_L007" = 10,
  "SIGAH4_E105_5_TET123_Chimera_TKO_L008" = 10
)
io$max_doublet_score$third_batch <- c(
  "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001" = 1000,
  "SIGAB6_E85_3_Dnmt3aWT_Dnmt3b_WT_L002" = 1000,
  "SIGAC6_E85_5_Dnmt3aKO_Dnmt3b_Het_L003" = 1000,
  "SIGAD6_E85_8_Dnmt3aHet_Dnmt3b_KO_L004" = 1000,
  "SIGAE4_E105_3_TET123_Chimera_Host_L005" = 1000,
  "SIGAF4_E105_3_TET123_Chimera_TKO_L006" = 1000,
  "SIGAG4_E105_5_TET123_Chimera_Host_L007" = 1000,
  "SIGAH4_E105_5_TET123_Chimera_TKO_L008" = 1000
)