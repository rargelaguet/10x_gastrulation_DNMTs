source(here::here("settings.R"))
source(here::here("utils.R"))
library(readxl)
library(cowplot)
library(org.Mm.eg.db)

io$celltype <- "ExE_ectoderm"#"Neural_crest"#"ExE_ectoderm"
dahlet_file <- "~/data/Dahlet_2020/sup_tables/41467_2020_16919_MOESM5_ESM.xlsx"

marker_genes <- fread(io$atlas.marker_genes) %>% 
  .[celltype == io$celltype] %>% 
  .[order(-rank(abs(logFC)))]

if (nrow(marker_genes) > 10) {
  markers <- marker_genes[1:10, gene]
} else {
  markers <- marker_genes[, gene]
}


d1_label <- "Dnmt1\nKO"

dahlet_exp <- read_xlsx(dahlet_file, sheet = 2) %>%
  setDT() %>%
  .[, c(1:7)] %>%
  setnames(colnames(.), gsub("#", "_", colnames(.))) %>%
  setnames("gene.id", "gene") %>%
  .[gene %in% markers] %>%
  melt(id.var = "gene", value.name = "exp", variable.name = "sample") %>%
  .[, type := strsplit(as.character(sample), "_") %>% map_chr(2)] %>%
  .[, type := factor(type, levels = c("WT", "D1KO"))] %>%
  .[, sample := gsub("fpkm_", "", sample)] %>%
  .[, sample := factor(sample, levels = c("WT_1", "WT_2", "WT_3", "D1KO_1", "D1KO_2", "D1KO_3"))] %>%
  .[, log2_exp := log2(exp+1)] %>%
  .[type == "WT", type2 := "WT"] %>%
  .[type == "D1KO", type2 := d1_label]

gene_aliases <- as.data.table(org.Mm.egALIAS2EG) %>% 
  merge(as.data.table(org.Mm.egENSEMBL), by = "gene_id") %>% 
  #merge(as.data.table(org.Mm.egSYMBOL), by = "gene_id") %>% 
  setnames(c("ensembl_id", "alias_symbol"), c("ens_id", "gene"))


dnmt1_dahlet <- read_xlsx(dahlet_file, sheet = 2) %>%
  setDT() %>%
  .[, .(gene = gene.id,
        dahlet_logFC = as.numeric(log2FoldChange),
        dahlet_padj = as.numeric(padj),
        gene.id.2
  )
  ] %>%
  .[, dahlet_sig := FALSE] %>%
  .[dahlet_padj < 0.001 & abs(dahlet_logFC) > 3, dahlet_sig := TRUE] %>% # cutoffs apparently used in the paper: padj < 0.001, logFC > 3
  .[dahlet_sig == T & dahlet_logFC > 1, dahlet_updown := "upregulated"] %>%
  .[dahlet_sig == T & dahlet_logFC < -1, dahlet_updown := "downregulated"] %>% 
  merge(gene_aliases, by = "gene", all.x = TRUE)





plot_bars <- function(genes){
  toplot <- dahlet_exp[gene %in% genes]
  to.plot.means <- toplot[,.(log2_exp=mean(log2_exp),sd=sd(log2_exp)), .(type, type2, gene)]
  
  max <- to.plot.means[, max(log2_exp)]
  
  p.vals <- dnmt1_dahlet[gene %in% genes] %>% 
    .[, c("group1", "group2", "y.position") := .("WT", d1_label, max * 1.2)] %>% 
    .[, p := paste("p =", signif(dahlet_padj, digits = 2))]
  
  palette <- c(
    "#FCF8CD",
    "#80B1D3"
    
  )
  
  ggplot(to.plot.means, aes(type2, log2_exp)) +
    geom_bar(stat = "identity", colour = "black", aes(fill = type2)) +
    geom_jitter(size=1, alpha=0.90, width=0.15, shape=21, data = toplot) +
    geom_errorbar(aes(ymin=log2_exp-sd, ymax=log2_exp+sd), width=0.25, alpha=0.75, size=0.5) +
    stat_pvalue_manual(data = p.vals, size = 2.5) +
    facet_wrap(~gene, nrow = 1) +
    theme_cowplot() +
    #guides(fill = FALSE, colour = FALSE) +
    #theme(strip.background =element_rect(fill="white")) +
    #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(
      strip.text = element_text(size=rel(0.85)),
      strip.background = element_blank(),
      axis.text.x = element_text(colour="black",size=rel(0.75)),
      axis.text.y = element_text(colour="black",size=rel(0.8)),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(colour="black",size=rel(1.0)),
      legend.position = "none"
    ) +
    labs(x = "Genotype", y = "log2 FPKM", fill = "") +
    scale_fill_manual(values=palette) 
}

(plot <- plot_bars(markers))
