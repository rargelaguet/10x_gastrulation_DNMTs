

plotting_fn2 <- function(dt, gene, classes, celltypes, max.expr=NULL) {
  
  dt.filt <- dt[repeat_class == gene & celltype.mapped%in%celltypes & class%in%classes]
  

  
  to.plot <- dt.filt[, .(
    cell = alias,
    expr,
    class = factor(class, levels = classes),
    celltype = factor(celltype.mapped, levels = celltypes)
  )]
  
  levels(to.plot$class) <- gsub("_", "\n", levels(to.plot$class))
  names(opts$classes.colors) <- gsub("_", "\n", names(opts$classes.colors))
  
  # For viz purposes
  if (!is.null(max.expr)) {
    to.plot[expr>=max.expr,expr:=max.expr]
  }
  to.plot.means <- to.plot[, .(expr = mean(expr), sd = sd(expr)), .(class, celltype)]
  
  
  ggplot(to.plot.means, aes(x=class, y=expr)) +
    geom_bar(stat = "identity", colour = "black", aes(fill = class)) +
    # geom_violin(scale = "width", alpha=0.8) +
    # geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.8) +
    # stat_summary(fun.data = give.n, geom = "text", size=3) +
    geom_jitter(size=1, shape=21, stroke=0.15, alpha=0.8, data=to.plot, width=0.1, aes(fill = class)) +
    geom_errorbar(aes(ymin=expr-sd, ymax=expr+sd), width=0.25, alpha=0.75, size=0.5) +
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

