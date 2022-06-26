plotting_fn <- function(sce, gene, classes, celltypes, max.expr=NULL) {
  
  sce.filt <- sce[gene,sce$celltype.mapped%in%celltypes & sce$class%in%classes]
  
  to.plot <- data.table(
    cell = colnames(sce.filt),
    expr = logcounts(sce.filt)[1,],
    class = factor(sce.filt$class, levels=classes),
    celltype = factor(sce.filt$celltype.mapped, levels=celltypes)
  )
  
  # For viz purposes
  if (!is.null(max.expr)) {
    to.plot[expr>=max.expr,expr:=max.expr]
  }
  to.plot.jitter <- to.plot[,.SD[sample.int(100)],by=c("class","celltype")]
  
  ggplot(to.plot, aes(x=class, y=expr, fill=class)) +
    geom_violin(scale = "width", alpha=0.8) +
    geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.8) +
    # stat_summary(fun.data = give.n, geom = "text", size=3) +
    geom_jitter(size=1, shape=21, stroke=0.15, alpha=0.4, data=to.plot.jitter, width=0.1) +
    scale_fill_manual(values=opts$classes.colors) +
    facet_wrap(~celltype, scales="fixed", nrow=1) +
    theme_classic() +
    labs(x="",y=sprintf("%s expression",gene)) +
    theme(
      strip.text = element_text(size=rel(0.85)),
      strip.background = element_blank(),
      axis.text.x = element_text(colour="black",size=rel(0.75)),
      axis.text.y = element_text(colour="black",size=rel(0.8)),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(colour="black",size=rel(1.0)),
      legend.position = "none"
    )
}