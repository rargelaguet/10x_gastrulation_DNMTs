source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")

to.merge <- c(
  "Erythroid3" = "Erythroid",
  "Erythroid2" = "Erythroid",
  "Erythroid1" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors",
  "Intermediate_mesoderm" = "Mixed_mesoderm",
  "Paraxial_mesoderm" = "Mixed_mesoderm",
  "Nascent_mesoderm" = "Mixed_mesoderm",
  "Pharyngeal_mesoderm" = "Mixed_mesoderm",
  "Visceral_endoderm" = "ExE_endoderm"
)

sample_metadata <- fread(io$metadata) %>% 
  .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped," ","_")] %>%
  .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,"/","_")] %>%
  .[,celltype.mapped2:=stringr::str_replace_all(celltype.mapped,to.merge)]

unique(sample_metadata$celltype.mapped)
unique(sample_metadata$celltype.mapped2)

fwrite(sample_metadata, io$metadata, sep="\t", col.names=T, row.names=F, na="NA", quote=F)
