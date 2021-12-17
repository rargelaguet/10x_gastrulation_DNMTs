io <- list()

io$min.counts.per.gene <- 50
io$k <- 30
io$min.mean <- 1e-3
io$npcs <- 50

io$atlas.RDS <- NULL # Choose between NULL and path to atlas SCE RDS
io$atlas.metadata <- NULL # Choose between NULL and path to atlas metadata txt
io$corrected.atlas.RDS <- NULL # Choose between NULL and path to atlas SCE RDS
io$corrected.atlas.metadata <- NULL  # Choose between NULL and path to atlas metadata txt
io$atlas_stages <- NULL # Choose between NULL and atlas stages
io$testing <- NULL # Choose between NULL and TRUE

io$subset_celltypes <- list(
    "Haematoendothelial" = c("Haematoendothelial progenitors", "Blood progenitors 1", "Blood progenitors 2", "Erythroid1", "Erythroid2", "Erythroid3", "Endothelium", "Cardiomyocytes"),
    "Endoderm" = c("ExE endoderm", "Visceral endoderm", "Gut", "Def. endoderm", "Notochord"),
    "EpiblastPSNeuro" = c("Epiblast", "Primitive Streak", "Rostral neurectoderm", "Spinal cord", "Surface ectoderm", "Nascent mesoderm", "NMP", "Neural crest", "Caudal neurectoderm", "Caudal epiblast", "Anterior Primitive Streak", "Forebrain/Midbrain/Hindbrain", "PGC"),
    "Mesoderm" = c("Nascent mesoderm", "Mesenchyme", "Mixed mesoderm", "ExE mesoderm, Intermediate mesoderm", "Pharyngeal mesoderm", "Paraxial mesoderm", "Somitic mesoderm", "Caudal Mesoderm", "Allantois")
    # Leaving out Parietal Endoderm and ExE ectoderm because they should already be well identified.
)

io$order <- list()
io$order[["first_batch"]] <- list()
io$order[["second_batch"]] <- list()
io$order[["third_batch"]] <- list()

io$order[["second_batch"]][["WT"]] <- c("ATLAS", "E85_Rep2_WT_Host_L005", "E85_Rep1_WT_Host_L003", "E75_WT_Host_L001")
io$order[["second_batch"]][["TET_TKO"]] <- c("ATLAS", "E85_Rep2_TET_TKO_L006", "E85_Rep1_TET_TKO_L004", "E75_TET_TKO_L002")
io$order[["second_batch"]][["DNMT3B_HET"]] <- c("ATLAS", "E125_DNMT3A_HET_A_L001", "E125_DNMT3A_HET_A_L003")
io$order[["second_batch"]][["DNMT3B_KO"]] <- c("ATLAS", "E125_DNMT3A_KO_E_L004", "E125_DNMT3A_KO_B_L002")

io$order[["third_batch"]][["TET_Host"]] <- c("ATLAS", "SIGAE4_E105_3_TET123_Chimera_Host_L005", "SIGAG4_E105_5_TET123_Chimera_Host_L007")
io$order[["third_batch"]][["TET_TKO"]] <- c("ATLAS", "SIGAF4_E105_3_TET123_Chimera_TKO_L006", "SIGAH4_E105_5_TET123_Chimera_TKO_L008")
io$order[["third_batch"]][["Dnmt3aKO_Dnmt3bWT"]] <- c("ATLAS", "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001")
io$order[["third_batch"]][["Dnmt3aWT_Dnmt3bWT"]] <- c("ATLAS", "SIGAB6_E85_3_Dnmt3aWT_Dnmt3b_WT_L002")
io$order[["third_batch"]][["Dnmt3aKO_Dnmt3bHet"]] <- c("ATLAS", "SIGAC6_E85_5_Dnmt3aKO_Dnmt3b_Het_L003")
io$order[["third_batch"]][["Dnmt3aHet_Dnmt3bKO"]] <- c("ATLAS", "SIGAD6_E85_8_Dnmt3aHet_Dnmt3b_KO_L004")