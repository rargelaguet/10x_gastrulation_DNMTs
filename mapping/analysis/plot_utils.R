colors <- c(
  "Ectoderm" = "steelblue",
  "Mesoderm" = "#CD3278",
  "Endoderm" = "#43CD80",
  "Blood" = "#C72228",

  "Epiblast" = "#635547",
  "Primitive Streak" = "#DABE99",
  "Caudal epiblast" = "#9e6762",

  "PGC" = "#FACB12",

  "Anterior Primitive Streak" = "#c19f70",
  "Notochord" = "#0F4A9C",
  "Def. endoderm" = "#F397C0",
  "Gut" = "#EF5A9D",

  "Nascent mesoderm" = "#C594BF",
  "Mixed mesoderm" = "#DFCDE4",
  "Intermediate mesoderm" = "#139992",
  "Caudal Mesoderm" = "#3F84AA",
  "Paraxial mesoderm" = "#8DB5CE",
  "Somitic mesoderm" = "#005579",
  "Pharyngeal mesoderm" = "#C9EBFB",
  "Cardiomyocytes" = "#B51D8D",
  "Allantois" = "#532C8A",
  "ExE mesoderm" = "#8870ad",
  "Mesenchyme" = "#cc7818",

  "Haematoendothelial progenitors" = "#FBBE92",
  "Endothelium" = "#ff891c",
  "Blood progenitors" = "#f9decf",
  "Blood progenitors 1" = "#f9decf",
  "Blood progenitors 2" = "#c9a997",
  "Erythroid1" = "#C72228",
  "Erythroid2" = "#f79083",
  "Erythroid" = "#f79083",
  "Erythroid3" = "#EF4E22",

  "NMP" = "#8EC792",

  "Rostral neurectoderm" = "#65A83E",
  "Caudal neurectoderm" = "#354E23",
  "Neural crest" = "#C3C388",
  "Forebrain/Midbrain/Hindbrain" = "#647a4f",
  "Spinal cord" = "#CDE088",

  "Surface ectoderm" = "#f7f79e",

  "Visceral endoderm" = "#F6BFCB",
  "ExE endoderm" = "#7F6874",
  "ExE ectoderm" = "#989898",
  "Parietal endoderm" = "#1A1A1A"
                    
)

plot.dimred <- function(plot_df) {
  ggplot(data=plot_df, mapping=aes(x=V1, y=V2)) +
    ggrastr::geom_point_rast(aes(size=mapped, alpha=mapped, colour=mapped)) +
    scale_size_manual(values = c("TRUE"=opts$size.mapped, "FALSE"=opts$size.nomapped)) +
    scale_alpha_manual(values = c("TRUE"=opts$alpha.mapped, "FALSE"=opts$alpha.nomapped)) +
    labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
    scale_colour_manual(values = c("TRUE"="red", "FALSE"="grey")) +
    guides(colour = guide_legend(override.aes = list(size=6))) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
}

plot.dimred.wtko <- function(plot_df) {
  ggplot(data=plot_df, mapping=aes(x=V1, y=V2)) +
    ggrastr::geom_point_rast(aes(size=mapped, alpha=mapped, colour=mapped)) +
    scale_size_manual(values = c("-10"=opts$size.mapped, "10"=opts$size.mapped, "0"=opts$size.nomapped)) +
    scale_alpha_manual(values = c("-10"=opts$alpha.mapped, "10"=opts$alpha.mapped, "0"=opts$alpha.nomapped)) +
    scale_colour_manual(values = c("-10"="red", "10"="blue", "0"="lightgrey")) +
    guides(colour = guide_legend(override.aes = list(size=6))) +
    theme_classic() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank()
    )
}
