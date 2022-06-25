library(data.table)
library(purrr)
library(ggplot2)
library(cowplot)

# this script makes boxplots of sites overlapping DEGs from our analysis

#metfile     <- "/Users/sclark/data/Dahlet_2020/wgbs/feature_level/Mmusculus_genes_BioMart.87.tsv.gz"
  #"/Users/sclark/data/Dahlet_2020/wgbs/feature_level/multiome_peaks.tsv.gz"
metfile <-  "/Users/sclark/data/Dahlet_2020/wgbs/feature_level/prom_2000_2000.tsv.gz"
 

genesfile   <- "/Users/sclark/data/Dahlet_2020/wgbs/annoations/Mmusculus_genes_BioMart.87.txt"

degenes <- "/Users/sclark/data/Dahlet_2020/10x_Dnmts_results/differential/Dnmt1_KO/"
fdr_threshold <- 0.01
logFC_threshold <- 1

extend_by <- 50000 # when finding overlapping genes extend each locus up and down by this many bp. Note thatif annotation contains ens_id, overlap is not performed.


sample_names <- data.table(
  sample = c(
    "GSM3752614_WGBS_WT_Rep1",
    "GSM3752615_WGBS_D1KO_Rep1",
    "GSM3752616_WGBS_DKO_Rep1",
    "GSM4558210_WGBS_WT_Rep2",
    "GSM4558211_WGBS_D1KO_Rep2",
    "GSM4558212_WGBS_DKO_Rep2"
  ),
  sample_name = c(
    "WT rep 1",
    "Dnmt1 KO rep 1",
    "Dnmt3 DKO rep 1",
    "WT rep 2",
    "Dnmt1 KO rep 2",
    "Dnmt3 DKO rep 2"
  )
) %>%
  .[, sample_name := factor(sample_name, levels = c(
    "WT rep 1",
    "Dnmt1 KO rep 1",
    "Dnmt3 DKO rep 1",
    "WT rep 2",
    "Dnmt1 KO rep 2",
    "Dnmt3 DKO rep 2"
  ))]



de <- dir(degenes, full=T, pattern = "txt.gz$") %>%
  map(~{
    type <- gsub(".txt.gz", "", basename(.x))
    dt <- fread(.x) %>%
      .[padj_fdr < fdr_threshold] %>%
      .[abs(logFC) >= logFC_threshold] %>% 
      .[, sign := sign(logFC)] %>%
      .[sign == -1, updown := "Upregulated in KO"] %>%
      .[sign == 1, updown := "Downregulated in KO"] %>%
      .[, type := type]
      #.[, .(ens_id, updown, type = type)]
  }) %>%
  rbindlist()

genes <- fread(genesfile) %>%
  setnames("symbol", "gene") %>%
  .[, c("start", "end") := .(start - extend_by, end + extend_by)] %>%
  setkey(chr, start, end)

met <- fread(metfile) %>%
  merge(sample_names, by = "sample") %>%
  setkey(chr, start, end)

if (met[1, id] %like% "ENSMUSG"){
  met[, ens_id := id]
  metgenes <- met
} else {
  metgenes <- foverlaps(genes, met, nomatch = 0L)
}



metde <- merge(metgenes, de, by = "ens_id", all.x = TRUE ,allow.cartesian = TRUE) %>%
  .[is.na(updown), updown := "Not DE"] %>%
  .[, updown := factor(updown, levels = c(
    "Not DE",
    "Downregulated in KO",
    "Upregulated in KO"
  ))] %>% 
  .[, de := "DEG"] %>% 
  .[updown == "Not DE", de := "Not DE"]

metde_allcelltypes <- metde[, .(ens_id, sample, anno, rate, N, sample_name, updown, de)] %>%
  unique()

# ggplot(metde[type %like% "Blood_progenitors"], aes(updown, rate, fill = updown, colour = updown)) +
#   geom_boxplot(alpha = 0.5, outlier.shape = NA) +
#   facet_wrap(~sample_name) +
#   theme_cowplot() +
#   guides(fill = FALSE, colour = FALSE) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(metde_allcelltypes[], aes(updown, rate, fill = updown, colour = updown)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  facet_wrap(~sample_name) +
  theme_cowplot() +
  guides(fill = FALSE, colour = FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Methylation at promoters of Dnmt1 KO DE genes")

ggplot(metde_allcelltypes[], aes(de, rate, fill = de, colour = de)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  facet_wrap(~sample_name) +
  theme_cowplot() +
  guides(fill = FALSE, colour = FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Methylation at promoters of Dnmt1 KO DE genes")

# ggplot(metde[type %like% "Rostral_neurectoderm"], aes(sample, rate, fill = sample, colour = sample)) +
#   geom_boxplot(alpha = 0.5, outlier.shape = NA) +
#   facet_wrap(~updown) +
#   theme_cowplot() +
#   guides(fill = FALSE, colour = FALSE) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))