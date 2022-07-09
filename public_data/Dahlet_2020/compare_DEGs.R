library(data.table)
library(purrr)
library(ggplot2)
library(readxl)
library(cowplot)
library(ggVennDiagram)
library(VennDiagram)
library(ggpubr)
library(org.Mm.eg.db)
library(ggrastr)
library(ggrepel)

# this script makes some comparisons between our DEGs and those reported by Dahlet et al 2020


dahlet_file <- "~/data/Dahlet_2020/sup_tables/41467_2020_16919_MOESM5_ESM.xlsx"
our_files <- "~/data/Dahlet_2020/10x_Dnmts_results/differential/"

plots_out_dir <- "/Users/sclark/Google Drive/My Drive/projects/babraham/dnmt_tet/plots/dahlet2020/exp/"

epi_genes <- c(
  "Pou5f1", "Utf1", "Slc7a3", "Fgf5", "Pim2"
)
exe_genes <- c(
  "Rhox5", "Krt8", "Apoe", "Ascl2", "Trap1a" , "Xlr3a"
)

hox_genes <- c(
  "Hoxc9", "Hoxc8", "Hoxb9" , "Hoxa9"
)

# load ensembl IDs to enable merging

gene_aliases <- as.data.table(org.Mm.egALIAS2EG) %>% 
  merge(as.data.table(org.Mm.egENSEMBL), by = "gene_id") %>% 
  #merge(as.data.table(org.Mm.egSYMBOL), by = "gene_id") %>% 
  setnames(c("ensembl_id", "alias_symbol"), c("ens_id", "gene"))


# from the paper: 'This analysis identified 414 upregulated and 68 downregulated genes in Dnmt1−/− embryos'
# also from the paper: 'differentially expressed genes were determined using DESeq2 v1.16.1 (fold change > 3, adjusted p value < 0.001)'
# they clearly didn't use these cutoffs. It looks like the 414 upregulated genes are the first 414 genes in the xlsx. Not sure about the downregulated genes. 
# its possible they are using some other cutoff (e.g. for total expression levels)
# but the top 414 genes are not all >3 logFC so who knows what's going on

# for this analysis I will just use the cutoffs they state and ignore the fact that we get different numbers of genes

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


dnmt1_dahlet[, table(dahlet_updown)]
dnmt1_dahlet[dahlet_sig == T][order(dahlet_updown, dahlet_padj)] 
dnmt1_dahlet[dahlet_sig == T][order(-dahlet_logFC)] 

# check all sig hits have an ensembl id
dnmt1_dahlet[dahlet_sig==T & is.na(ens_id)]

dnmt1_dahlet[dahlet_sig==T & is.na(ens_id), .N]
dnmt1_dahlet[dahlet_sig==T, .N]

# 6 out of 306 missing... 

dnmt1_ours <- dir(paste0(our_files, "/Dnmt1_KO"), f = T, pattern = ".txt.gz$") %>% 
  set_names(gsub(".txt.gz$", "", basename(.))) %>% 
  map(fread) %>% 
  map2(names(.), ~.x[, type := gsub("_WT.*", "", .y)]) %>% 
  rbindlist() %>% 
  .[, type2 := gsub("Forebrain_Midbrain_Hindbrain","Forebrain / Midbrain / Hindbrain", type) %>% gsub("_", " ", .)] %>% 
  .[, sig := ifelse(padj_fdr < 0.01 & abs(logFC) >=1, TRUE, FALSE)] # cutoffs we used in our anlaysis: padj < 1%, logFC >= 1

setnames(dnmt1_ours, "gene", "gene_ours")

overlap <- merge(dnmt1_dahlet, dnmt1_ours, all = TRUE, by = "ens_id")


overlap[dahlet_sig == T]

overlap[dahlet_sig == T & dahlet_updown == "upregulated", .(mean(logFC, na.rm = T), mean(padj_fdr, na.rm = T)), type]
overlap[dahlet_sig == T & dahlet_updown == "downregulated", .(mean(logFC, na.rm = T), mean(padj_fdr, na.rm = T)), type]

overlap[, dot_alpha := 0.5]
overlap[abs(logFC)>2.5 | abs(dahlet_logFC)>2.5, dot_alpha := 2]

overlap[, show_gene := FALSE]
overlap[gene %in% c(exe_genes, hox_genes, epi_genes) & sig == TRUE, show_gene := TRUE]

# reduce the number of dots to plot

toplot <- rbind(
  overlap[abs(logFC)<1][sample.int(.N/2)],
  overlap[abs(logFC)>=1]
)


types <- overlap[, unique(type)]
types_sub <- types[1:16]
types_sub=types

types_sub <- c(
  # "Gut",
  # "ExE_ectoderm" ,
  # "Rostral_neurectoderm",
  # "Haematoendothelial_progenitors",
  "Surface_ectoderm"
)

(scatterplot <- ggplot(toplot[type %in% types_sub], aes(dahlet_logFC, logFC)) +
  geom_point(size = 0.5, aes(alpha = dot_alpha)) +
  geom_smooth(method = "lm", colour = "black", size = 0.5) +
  geom_text_repel(data = toplot[type %in% types_sub][show_gene == TRUE], aes(label = gene)) + 
  theme_cowplot() +
  facet_wrap(~type2, ncol = 1) +
  theme(strip.background =element_rect(fill="white")) +
  guides(alpha = FALSE) +
  labs(x = "Dahlet et al 2020 (log2 fold change)", y = "This Study (log2 fold change)") +
  stat_cor(method = "pearson")
  )

(scatterplot <- ggplot(toplot[], aes(dahlet_logFC, logFC)) +
    geom_point(size = 0.5, aes(alpha = dot_alpha)) +
    geom_smooth(method = "lm", colour = "black", size = 0.5) +
    geom_text_repel(data = toplot[][show_gene == TRUE], aes(label = gene)) + 
    theme_cowplot() +
    facet_wrap(~type2, ncol = 5) +
    theme(strip.background =element_rect(fill="white")) +
    guides(alpha = FALSE) +
    labs(x = "Dahlet et al 2020 (log2 fold change)", y = "This Study (log2 fold change)") +
    stat_cor(method = "pearson")
)

# (scatterplot <- ggplot(overlap[], aes(dahlet_logFC, logFC)) +
#     geom_point(alpha=0.1, size = 0.5) +
#     geom_smooth(method = "lm", colour = "black", size = 0.5) +
#     theme_cowplot() +
#     facet_wrap(~type2, ncol = 1) +
#     theme(strip.background =element_rect(fill="white")) +
#     labs(x = "Dahlet et al 2020", y = "This Study") +
#     stat_cor(method = "pearson")
# )



(scatterplot <- rasterize(scatterplot, layers='Point', dpi=150))

outfile <- "/Users/sclark/Google Drive/My Drive/projects/babraham/dnmt_tet/plots/dahlet2020/exp/scatter.pdf"

save_plot(outfile, scatterplot, base_height = 25, base_width = 25)

# same plots as individual pdfs

scatterplots <- map(types[!is.na(types)], ~{
  p <- ggplot(toplot[type %in% .x], aes(dahlet_logFC, logFC)) +
     geom_point(size = 0.5, aes(alpha = dot_alpha)) +
     geom_smooth(method = "lm", colour = "black", size = 0.5) +
     geom_text_repel(data = toplot[type %in% .x][show_gene == TRUE], aes(label = gene)) + 
     theme_cowplot() +
     facet_wrap(~type2, ncol = 1) +
     theme(strip.background =element_rect(fill="white")) +
     guides(alpha = FALSE) +
     labs(x = "Dahlet et al 2020 (log2 fold change)", y = "This Study (log2 fold change)") +
     stat_cor(method = "pearson")
  
  rasterize(p, layers='Point', dpi=150)
})
 outfiles <- paste0(
   "/Users/sclark/Google Drive/My Drive/projects/babraham/dnmt_tet/plots/dahlet2020/exp/scatter/",
   types[!is.na(types)],
   ".pdf"
 )
 dir.create(dirname(outfiles[1]))
 
 walk2(outfiles, scatterplots, save_plot)

# what is the overlap in sig hits?

# some gene ids do not match between the two datasets

# (missing_in_ours <- overlap[dahlet_sig==T & is.na(ens_id), unique(gene)])
# missing_in_dahlet <- overlap[is.na(gene.id.2), unique(gene)]
# overlap[dahlet_sig == T]
# overlap[gene != gene.id.2]
# 
# overlap[is.na(gene.id.2)]



venn <- list(
  `Dahlet et al` = overlap[dahlet_sig == T & !is.na(ens_id), ens_id],
  `This study`   = overlap[sig == T & !is.na(gene.id.2), unique(ens_id)]
)


(vennplot <- ggVennDiagram(
  venn, 
  label = "count", 
  label_alpha = 0,
  edge_size = 0.2
  ) +
  guides(fill=F) +
  scale_fill_gradient(low="white",high = "white")
)

outfile <- paste0(plots_out_dir, "/venn.pdf")
save_plot(outfile, vennplot)



# find examples of genes DE in our analysis but not in the bulk
overlap[sig ==T & dahlet_sig == F][order(-rank(abs(logFC)))]




# subset to genes we mention in the text



overlap[sig ==T & dahlet_sig == F][gene %in% c(epi_genes, exe_genes, hox_genes)][,table(sign(logFC))]

overlap[sig ==T & dahlet_sig == F][gene %in% c(epi_genes, exe_genes, hox_genes)][order(-rank(abs(logFC)))]
overlap[sig ==T & dahlet_sig == F][gene %in% c(epi_genes, exe_genes, hox_genes)][order(-rank(abs(logFC)))][1:30]

hits <- overlap[sig ==T & dahlet_sig == F][gene %in% c(epi_genes, exe_genes, hox_genes)] 

hits[, .N, gene][order(-rank(N))]

#  top hits are Rhox5 and Trap1a so lets plot these

genes_toplot <- c(
  "Rhox5",
  "Trap1a"#,
  #"Xlr3a"
)

dahlet_exp <- read_xlsx(dahlet_file, sheet = 2) %>%
  setDT() %>% 
  .[, c(1:7)] %>% 
  setnames(colnames(.), gsub("#", "_", colnames(.))) %>% 
  setnames("gene.id", "gene") %>% 
  melt(id.var = "gene", value.name = "exp", variable.name = "sample") %>% 
  .[, type := strsplit(as.character(sample), "_") %>% map_chr(2)] %>% 
  .[, type := factor(type, levels = c("WT", "D1KO"))] %>% 
  .[, sample := gsub("fpkm_", "", sample)] %>% 
  .[, sample := factor(sample, levels = c("WT_1", "WT_2", "WT_3", "D1KO_1", "D1KO_2", "D1KO_3"))] %>% 
  .[, log2_exp := log2(exp+1)]

ggplot(dahlet_exp[gene %in% genes_toplot], aes(sample, log2_exp, colour = type, fill = type)) +
  geom_bar(stat = "identity") +
  facet_wrap(~gene) +
  theme_cowplot() +
  guides(fill = FALSE, colour = FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# for a sup fig, lets plot all the same genes we plot in fig2:

hox_genes <- c(
  "Hoxc9", 
  "Hoxc8", 
  "Hoxb9" , 
  "Hoxa9"
)

# hox_genes <- c(
#   "Hoxc9" ,
#   "Hoxa9"
# )


pluri_genes <- c(
  "Utf1",
  "Slc7a3",
  "Pou5f1",
  "Pim2",
  "Nanog",
  "Gng3",
  "Fgf3"
)

# pluri_genes <- c(
#   
#   "Pou5f1",
#   "Pim2"
#   
# )


exe_genes <-c(
  "Xlr3a",
  "Trap1a",
  "Tex19.1",
  "Rhox9",
  "Rhox5",
  "Fmr1nb",
  "Ascl2",
  "Apoe"
)

# exe_genes <-c(
#   "Trap1a",
#   "Rhox9"
#   
# )

plot_bars_by_sample <- function(genes){
  
  toplot <- dahlet_exp[gene %in% genes]
  
  palette <- c(
    "#FCF8CD",
    "#80B1D3"
  )
  
  ggplot(toplot, aes(sample, log2_exp, colour = type, fill = type)) +
    geom_bar(stat = "identity", colour = "black") +
    facet_wrap(~gene) +
    theme_cowplot() +
    #guides(fill = FALSE, colour = FALSE) +
    theme(strip.background =element_rect(fill="white")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Sample", y = "log2 FPKM", fill = "") +
    scale_fill_manual(values=palette)
}



plot_bars <- function(genes){
  toplot <- dahlet_exp[gene %in% genes]
  to.plot.means <- toplot[,.(log2_exp=mean(log2_exp),sd=sd(log2_exp)), .(type, gene)]
  
  max <- to.plot.means[, max(log2_exp)]
  
  p.vals <- dnmt1_dahlet[gene %in% genes] %>% 
    .[, c("group1", "group2", "y.position") := .("WT", "D1KO", max * 1.2)] %>% 
    .[, p := paste("p =", format(dahlet_padj, digits = 3))]
  
  palette <- c(
    "#FCF8CD",
    "#80B1D3"
    
  )
  
  ggplot(to.plot.means, aes(type, log2_exp)) +
    geom_bar(stat = "identity", colour = "black", aes(fill = type)) +
    geom_jitter(size=1, alpha=0.90, width=0.15, shape=21, data = toplot) +
    geom_errorbar(aes(ymin=log2_exp-sd, ymax=log2_exp+sd), width=0.25, alpha=0.75, size=0.5) +
    stat_pvalue_manual(data = p.vals, size = 2.5) +
    facet_wrap(~gene, nrow = 1) +
    theme_cowplot() +
    guides(fill = FALSE, colour = FALSE) +
    theme(strip.background =element_rect(fill="white")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Genotype", y = "log2 FPKM", fill = "") +
    scale_fill_manual(values=palette) 
}

(hoxplot <- plot_bars(hox_genes))
(pluriplot <- plot_bars(pluri_genes))
(exeplot <- plot_bars(exe_genes))

plots <- map(list(hox_genes, pluri_genes, exe_genes), plot_bars)







outfiles <- paste0(
  "/Users/sclark/Google Drive/My Drive/projects/babraham/dnmt_tet/plots/dahlet2020/exp/",
  c("hox.pdf", "pluri.pdf", "exe.pdf")
  
)
dirname(outfiles)[1] %>% dir.create(recursive = T)

walk2(outfiles, plots, save_plot)



# results of Dahlet's DEG analyses on these genes:

dnmt1_dahlet[gene %in% c(exe_genes, pluri_genes, hox_genes)]




