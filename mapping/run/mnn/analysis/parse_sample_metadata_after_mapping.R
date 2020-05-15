#########
## I/O ##
#########

source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")

io$mapping <- "/Users/ricard/data/10x_gastrulation_DNMTs/results/second_batch/mapping/mapping10x_mnn.rds"

###############
## Load data ##
###############

# Load mapping results
foo <- readRDS(io$mapping)$mapping %>% .[,c("cell","celltype.mapped")] %>%
  as.data.table %>% .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,"_"," ")]

bar <- fread("/Users/ricard/data/10x_gastrulation_DNMTs/results/second_batch/mapping/leah/mapping_mnn_WT.txt.gz") %>%
  .[cell %in% mapping.dt$cell] %>% .[,c("cell","celltype.mapped")] %>%
  .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,"/"," ")]

foobar <- merge(foo,bar,by="cell")

mean(foobar$celltype.mapped.x == foobar$celltype.mapped.y)

disagreement <- foobar[!celltype.mapped.x == celltype.mapped.y]

fwrite(disagreement, "/Users/ricard/data/10x_gastrulation_DNMTs/results/second_batch/mapping/disagreement_leah_vs_ricard.txt.gz")

#################
## Save output ##
#################

fwrite(sample.metadata, file=io$sample_metadata, sep="\t", col.names=T, row.names=F, na="NA", quote=F)


