celltype_colours = c(
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
 "Blood progenitors 1" = "#f9decf",
 "Blood progenitors 2" = "#c9a997",
 "Erythroid1" = "#C72228",
 "Erythroid2" = "#f79083",
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

barplot.pub <- function(to.plot, x="lineage10x_2", colors=NULL, xlim.max=NULL) {
  p <- ggplot(to.plot, aes_string(x=x, y="N")) +
    scale_x_discrete(drop=FALSE) + 
    # coord_flip() +
    labs(y="Number of cells") +
    theme_classic() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=rel(1.3)),
      axis.title.x = element_text(color="black", size=rel(1.1)),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1.3), color="black"),
      axis.text.x = element_text(size=rel(1.1), color="black")
    )
    
    if (is.null(colors)) {
        p <- p + geom_bar(stat="identity", color="black")
    } else {
        p <- p + geom_bar(aes_string(fill=x), stat="identity", color="black") + 
            scale_fill_manual(values=colors, drop=F)
    }

    if (!is.null(xlim.max)) {
      p <- p + coord_flip(ylim=c(0,xlim.max))
    } else {
      p <- p + coord_flip()
    }
  
    return(p)
}

stage_colours <- c("E6.5" = "#D53E4F",
                   "E6.75" = "#F46D43",
                   "E7.0" = "#FDAE61",
                   "E7.25" = "#FEE08B",
                   "E7.5" = "#FFFFBF",
                   "E7.75" = "#E6F598",
                   "E8.0" = "#ABDDA4",
                   "E8.25" = "#66C2A5",
                   "E8.5" = "#3288BD",
                   "mixed_gastrulation" = "#A9A9A9")

sample_colours <- c("ESC" = "#635547",
                    "nomap" = "gray",
                    "EB_d3" = "#cfd9e4",
                    "EB_d4" = "#88a1bc",
                    "EB_d5" = "#406993",
                    "EB_d6" = "#2e4660",
                    "gastr.d3" = "#88bc88",
                    "gastr.d4" = "#409341",
                    "MPP_d2" = "#bc88bb",
                    "MPP_d3" = "#934092",
                    "MPA_d2" = "#c28873",
                    "MPA_d3" = "#ae664c"
                    )

io <- list()

# Dot size
io$size.mapped <- 0.6
io$size.nomapped <- 0.1
# for atlas only
io$dot_size = 0.6

# Transparency
io$alpha.mapped <- 1.0
io$alpha.nomapped <- 0.5
# for atlas only
io$dot_alpha = 1

io$plot.type <- "pdf" # NULL or one of "eps", "ps", "tex" (pictex), "pdf", "jpeg", "tiff", "png", "bmp", "svg"