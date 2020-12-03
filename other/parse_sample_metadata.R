source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")

sample_metadata <- fread(io$metadata)

sample_metadata[,class:=stringr::str_replace_all(batch,opts$batch.to.class)]
unique(sample_metadata$class)

sample_metadata[,stage:="E8.5"] %>% .[grep("E12.5|E125|E12_5",batch),stage:="E12.5"]
unique(sample_metadata$stage)

fwrite(sample_metadata, io$metadata, sep="\t", quote=F, na="NA")
