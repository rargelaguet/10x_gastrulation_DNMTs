source(here::here("settings.R"))

#####################
## Define settings ##
#####################

# I/O
io$outdir <- paste0(io$basedir,"/results_new/qc/figure"); dir.create(io$outdir, showWarnings = F)

# Options
opts$min_nFeature_RNA <- 1000
opts$max_nFeature_RNA <- 10000
opts$mit_percent_RNA <- 15
opts$rib_percent_RNA <- 40

###################
## Load metadata ##
###################

metadata <- fread(io$metadata) %>% 
    .[class%in%opts$classes] %>% .[,class:=factor(class,levels=rev(opts$classes))] %>%
    .[,pass_rnaQC:=nFeature_RNA<=opts$max_nFeature_RNA & nFeature_RNA>=opts$min_nFeature_RNA & mit_percent_RNA<opts$mit_percent_RNA & rib_percent_RNA<opts$rib_percent_RNA]

###########################################
## Boxplots of QC metrics after QC calls ##
###########################################

qc.dt <- metadata %>% .[pass_rnaQC==TRUE] %>%
    melt(id.vars=c("alias","class","cell","stage"), measure.vars=c("nFeature_RNA","mit_percent_RNA","rib_percent_RNA"))

tmp <- data.table(
    variable = c("nFeature_RNA", "nFeature_RNA", "mit_percent_RNA", "rib_percent_RNA"),
    value = c(opts$min_nFeature_RNA, opts$max_nFeature_RNA, opts$mit_percent_RNA, opts$rib_percent_RNA)
)

# classes.to.plot <- c("Dnmt3a_KO", "WT", "Dnmt3ab_KO", "Dnmt3b_KO", "Dnmt1_KO")

# classes.to.plot <- unique(to.plot$class)
# for (i in classes.to.plot) {
#     
#     to.plot.jitter <- qc.dt[class==i] %>% .[sample.int(n=nrow(.), size=nrow(.)/3)]
#     
#     p <- ggplot(qc.dt[class==i], aes_string(x="alias", y="value")) +
#         ggrastr::geom_jitter_rast(alpha=0.5, width=0.15, size=0.05, data=to.plot.jitter) +
#         geom_boxplot(fill="gray70", outlier.shape=NA, coef=1, alpha=0.9) +
#         geom_hline(aes(yintercept=value), linetype="dashed", data=tmp) +
#         facet_wrap(~variable, scales="free_y", nrow=1, labeller = as_labeller(facet.labels)) +
#         # scale_fill_manual(values=opts$stage.colors) +
#         labs(x="", y="") +
#         theme_classic() +
#         theme(
#             axis.text.y = element_text(colour="black",size=rel(0.9)),
#             axis.text.x = element_text(colour="black",size=rel(0.55), angle=20, hjust=1, vjust=1),
#             axis.title.x = element_blank()
#         )
#     
#     pdf(file.path(io$outdir,sprintf("%s_qc_metrics_boxplot.pdf",i)), width=8, height=4)
#     # pdf(sprintf("%s/qc_metrics_boxplot.pdf",io$outdir))
#     print(p)
#     dev.off()
# }

facet.labels <- c("nFeature_RNA" = "Num. of genes", "mit_percent_RNA" = "Mitochondrial %", "rib_percent_RNA" = "Ribosomal %")

to.plot <- qc.dt %>% setorder(-class) %>% .[,alias:=factor(alias,levels=unique(.[["alias"]]))]
to.plot.jitter <- to.plot %>% .[sample.int(n=nrow(.), size=nrow(.)/10)]

p <- ggplot(to.plot, aes_string(x="alias", y="value", fill="class")) +
    ggrastr::geom_jitter_rast(aes(color=class), alpha=0.5, width=0.15, size=0.05, data=to.plot.jitter) +
    geom_boxplot(outlier.shape=NA, coef=1, alpha=0.9) +
    geom_hline(aes(yintercept=value), linetype="dashed", data=tmp) +
    facet_wrap(~variable, scales="free_y", labeller = as_labeller(facet.labels), nrow=1) +
    # scale_fill_brewer(palette="Dark2") +
    # scale_color_brewer(palette="Dark2") +
    scale_fill_manual(values=opts$classes.colors) +
    scale_color_manual(values=opts$classes.colors) +
    guides(x = guide_axis(angle = 90)) +
    labs(x="", y="") +
    theme_classic() +
    theme(
        legend.title = element_blank(),
        legend.position = "top",
        axis.text.y = element_text(colour="black",size=rel(0.9)),
        # axis.text.x = element_text(colour="black",size=rel(0.55), angle=20, hjust=1, vjust=1),
        axis.text.x = element_text(colour="black",size=rel(0.45)),
        axis.title.x = element_blank()
    )

pdf(file.path(io$outdir,"qc_metrics_boxplot.pdf"), width=13, height=5)
print(p)
dev.off()

#########################################################
## Plot fraction of cells that pass QC for each sample ##
#########################################################

to.plot <- metadata %>%
    .[,mean(pass_rnaQC,na.rm=T),by=c("alias","class","stage")]

p <- ggbarplot(to.plot, x="alias", y="V1", fill="class") +
    # scale_fill_manual(values=opts$stage.colors) +
    labs(x="", y="Fraction of cells that pass QC (RNA)") +
    # scale_fill_brewer(palette="Dark2") +
    scale_fill_manual(values=opts$classes.colors) +
    facet_wrap(~class, nrow=1, scales="free_x") +
    theme(
        legend.position = "none",
        axis.text.y = element_text(colour="black",size=rel(0.8)),
        axis.text.x = element_text(colour="black",size=rel(0.50), angle=20, hjust=1, vjust=1),
    )

pdf(file.path(io$outdir,"qc_metrics_barplot.pdf"), width=8, height=6)
print(p)
dev.off()

