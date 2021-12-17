#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(optparse))
#suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(cowplot))
#suppressPackageStartupMessages(library(gridExtra))
#suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(umap))
suppressPackageStartupMessages(library(ggpubr))

###################
## Parse options ##
###################
 
option_list = list(
    make_option(c("-R", "--input.rds"), type="character", default=NULL, 
                help="path to the mapping output", metavar="character"),
    make_option(c("-m", "--input.meta"), type="character", default=NULL, 
                help="path to the mapping metadata", metavar="character"),
    make_option(c("-q", "--query.metadata"), type="character", default=NULL, 
                help="path to the query metadata", metavar="character"),
    make_option(c("-a", "--atlas.metadata"), type="character", default=NULL, 
                help="path to the atlas metadata", metavar="character"),
    make_option(c("-S", "--settings"), type="character", default=NULL, 
                help="path to settings file", metavar="character"),
    make_option(c("-O", "--output.base"), type="character", default=NULL, 
                help="path to output RDS file", metavar="character"),
    make_option(c("-p", "--plot.type"), type="character", default="pdf", 
                help="what device should be used in ggsave", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opts = parse_args(opt_parser);

if (is.null(opts$input.rds)){
    print_help(opt_parser)
    stop("The path to a mapping output must be supplied.n", call.=FALSE)
} else if (is.null(opts$input.meta)){
    print_help(opt_parser)
    stop("The path to a mapping metadata must be supplied.n", call.=FALSE)
} else if (is.null(opts$query.metadata)){
    print_help(opt_parser)
    stop("The path to the query metadata must be supplied.n", call.=FALSE)
} else if (is.null(opts$atlas.metadata)){
    print_help(opt_parser)
    stop("The path to the atlas.metadata metadata must be supplied.n", call.=FALSE)
} else if (is.null(opts$settings)) {
    print_help(opt_parser)
    stop("A settings file must be supplied.n", call.=FALSE)
} else if (is.null(opts$output.base)) {
    print_help(opt_parser)
    stop("A path to an RDS output file must be supplied.n", call.=FALSE)
}

source(opts$settings)

if (!is.null(io$plot.type)) opts$plot.type <- io$plot.type

message("Loading input...")

# Load atlas metadata
if (file_ext(opts$atlas.metadata) == "rds") {
    meta_atlas <- readRDS(paste0(opts$atlas.metadata)) %>% as.data.table
} else {
    meta_atlas <- fread(opts$atlas.metadata, sep="\t", na="NA", quote=F) %>% as.data.table
}

# Load query metadata
meta_mapping <- fread(opts$input.meta)
meta_mapping <- meta_mapping[which(!is.na(meta_mapping$celltype.mapped)),]

# Load mapping output
merge_output <- readRDS(paste0(opts$input.rds))


find_top2_celltypes <- function(row) {
    tbl <- table(unlist(row))
    top2 <- sort(-tbl)[c(1,2)]
    return(names(top2))
}

get_proportion <- function(row) {
    foo <- row["N"]/inst[which((inst["V1"]==row["V1"]) && (inst["batch"]==row["batch"])),N]
    print(foo)
    return(row)
}

celltypes.mapped <- as.data.table(merge_output$celltypes.mapped)
top2_celltypes <- as.data.table(t(apply(celltypes.mapped,1,find_top2_celltypes )))
top2_celltypes$cell <- names(merge_output$celltypes.mapped[[1]])
top2_celltypes$batch <- meta_mapping[match(names(merge_output$celltypes.mapped[[1]]),meta_mapping$cell),batch]
toplot <- top2_celltypes[, .(N = length(unique(cell))), by = c("V1","V2","batch")]
inst <- top2_celltypes[,.N, by=c("V1","batch")]
for (i in c(1:dim(toplot)[1])) {
    row <- toplot[i,]
    toplot$prop[i] <- row$N/inst[intersect(which(inst$V1 == row$V1), which(inst$batch==row$batch)),N]
}

hm_prop <- ggplot(toplot, aes(V1, V2, fill= prop)) + 
  facet_wrap(~batch, ncol=2) +
  geom_tile() +
  scale_fill_gradient(low="gray90", high="#2171b5") +
  theme_bw() + 
  theme(text = element_text(size=12),
        axis.text.x=element_text(angle=90, hjust=1)) +
  labs(x = "Top Celltype", y = "Second-best Celltype", fill = "Proportion")
ggsave(paste0(opts$output.base, "_celltypes_heatmap_prop.", opts$plot.type), plot = hm_prop, device=opts$plot.type, width = 15, height = 15, units = "in")

hm_N <- ggplot(toplot, aes(V1, V2, fill= N)) + 
  facet_wrap(~batch, ncol=2) +
  geom_tile() +
  scale_fill_gradient(low="gray90", high="#2171b5") +
  theme_bw() + 
  theme(text = element_text(size=12),
        axis.text.x=element_text(angle=90, hjust=1)) +
  labs(x = "Top Celltype", y = "Second-best Celltype", fill = "Occurrences")
ggsave(paste0(opts$output.base, "_celltypes_heatmap_N.", opts$plot.type), plot = hm_N, device=opts$plot.type, width = 15, height = 15, units = "in")

# get joint corrected matrix and rename cells according to origin
message("Joining datasets...")

correct_map <- merge_output$correct_map
correct_atlas <- merge_output$correct_atlas
rownames(correct_map) <- lapply(rownames(correct_map), function(x) paste(paste(strsplit(x, "_")[[1]][1:2], collapse = "_"), "query", sep="_"))
rownames(correct_atlas) <- lapply(rownames(correct_atlas), function(x) paste(paste(strsplit(x, "_")[[1]][1:2], collapse = "_"), "atlas", sep="_"))

meta_atlas$cell <- lapply(meta_atlas$cell, function(x) paste(x, "atlas", sep="_"))
meta_atlas$MAP <- 'ATLAS'
meta_mapping$cell <- lapply(meta_mapping$cell, function(x) paste(x, "query", sep="_"))
meta_mapping$MAP <- 'QUERY'
                                  
correct_joined <- rbind(correct_atlas,correct_map)
meta_joined <- rbindlist(list(meta_atlas, meta_mapping), fill=TRUE)

message("Running UMAP...")
# calculate joint UMAP
UMAP_joined <- umap(correct_joined, method='umap-learn')


message("Creating dataframes to be plotted...")
# create data frames to be plotted
tmp <- data.frame(UMAP_joined$layout)
tmp$cell <- rownames(tmp)
plot_df_joined <- tmp %>% merge(meta_joined, by=c("cell"))
plot_df_joined$celltype[(which(is.na(plot_df_joined$celltype)))] <- plot_df_joined$celltype.mapped[(which(is.na(plot_df_joined$celltype)))]

plot_df_joined$sample[(which(is.na(plot_df_joined$sample)))] <- plot_df_joined$batch[(which(is.na(plot_df_joined$sample)))]

plot_df_joined$cellstage <- plot_df_joined$stage
plot_df_joined$cellstage[which(is.na(plot_df_joined$cellstage))] <- plot_df_joined$stage.mapped[which(is.na(plot_df_joined$cellstage))]
plot_df_joined$stage[which(is.na(plot_df_joined$stage))] <- plot_df_joined$batch[which(is.na(plot_df_joined$stage))]
plot_df_joined <- plot_df_joined[order(plot_df_joined$MAP),]

plot_df_atlas <- plot_df_joined[which(plot_df_joined$MAP == 'ATLAS'),]
plot_df_query <- plot_df_joined[which(plot_df_joined$MAP == 'QUERY'),]

message("Plotting...")
joint_stage <- ggplot(data=plot_df_joined, mapping = aes(x=X1, y=X2, colour=stage)) +
  geom_point(size=io$dot_size, alpha=io$dot_alpha) +
  # ggrastr::geom_point_rast(size=io$dot_size, alpha=io$dot_alpha) +
  labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  theme_classic()
ggsave(paste0(opts$output.base, "_joint_stage.", opts$plot.type), plot = joint_stage, device=opts$plot.type, width = 7, height = 5, units = "in")

query_stage <- ggplot(data=plot_df_query, mapping = aes(x=X1, y=X2, colour=stage)) +
  geom_point(size=io$dot_size, alpha=io$dot_alpha) +
  # ggrastr::geom_point_rast(size=io$dot_size, alpha=io$dot_alpha) +
  labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  theme_classic()
ggsave(paste0(opts$output.base, "_query_stage.", opts$plot.type), plot = query_stage, device=opts$plot.type, width = 7, height = 5, units = "in")

atlas_stage <- ggplot(data=plot_df_atlas, mapping = aes(x=X1, y=X2, colour=stage)) +
  geom_point(size=io$dot_size, alpha=io$dot_alpha) +
  # ggrastr::geom_point_rast(size=opts$dot_size, alpha=opts$dot_alpha) +
  labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
  scale_colour_manual(values = stage_colours, name = "Stage") +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  theme_classic()
ggsave(paste0(opts$output.base, "_atlas_stage.", opts$plot.type), plot = atlas_stage, device=opts$plot.type, width = 7, height = 5, units = "in")

joint_sample <- ggplot(data=plot_df_joined, mapping = aes(x=X1, y=X2, colour=sample)) +
  geom_point(size=io$dot_size, alpha=io$dot_alpha) +
  # ggrastr::geom_point_rast(size=io$dot_size, alpha=io$dot_alpha) +
  labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  theme_classic()
ggsave(paste0(opts$output.base, "_joint_sample.", opts$plot.type), plot = joint_sample, device=opts$plot.type, width = 7, height = 5, units = "in")

query_sample <- ggplot(data=plot_df_query, mapping = aes(x=X1, y=X2, colour=batch)) +
  geom_point(size=io$dot_size, alpha=io$dot_alpha) +
  # ggrastr::geom_point_rast(size=io$dot_size, alpha=io$dot_alpha) +
  labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  theme_classic()
ggsave(paste0(opts$output.base, "_query_sample.", opts$plot.type), plot = query_sample, device=opts$plot.type, width = 7, height = 5, units = "in")

atlas_sample <- ggplot(data=plot_df_atlas, mapping = aes(x=X1, y=X2, colour=sample)) +
  geom_point(size=io$dot_size, alpha=io$dot_alpha) +
  # ggrastr::geom_point_rast(size=opts$dot_size, alpha=opts$dot_alpha) +
  labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  theme_classic()
ggsave(paste0(opts$output.base, "_atlas_sample.", opts$plot.type), plot = atlas_sample, device=opts$plot.type, width = 7, height = 5, units = "in")

atlas_celltype <- ggplot(data=plot_df_atlas, mapping = aes(x=X1, y=X2, colour=celltype)) +
  geom_point(size=io$dot_size, alpha=io$dot_alpha) +
  # ggrastr::geom_point_rast(size=opts$dot_size, alpha=opts$dot_alpha) +
  labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
  scale_colour_manual(values = celltype_colours, name = "Celltype") +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  theme_classic()
ggsave(paste0(opts$output.base, "_atlas_celltype.", opts$plot.type), plot = atlas_celltype, device=opts$plot.type, width = 10, height = 5, units = "in")


query_celltype <- ggplot(data=plot_df_query, mapping = aes(x=X1, y=X2, colour=celltype.mapped)) +
  geom_point(size=io$dot_size, alpha=io$dot_alpha) +
  # ggrastr::geom_point_rast(size=opts$dot_size, alpha=opts$dot_alpha) +
  labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
  scale_colour_manual(values = celltype_colours, name = "Celltype") +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  theme_classic()
ggsave(paste0(opts$output.base, "_query_celltype.", opts$plot.type), plot = query_celltype, device=opts$plot.type, width = 10, height = 5, units = "in")


joint_celltype <- ggplot(data=plot_df_joined, mapping = aes(x=X1, y=X2, colour=celltype)) +
  geom_point(size=io$dot_size, alpha=io$dot_alpha) +
  # ggrastr::geom_point_rast(size=opts$dot_size, alpha=opts$dot_alpha) +
  labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
  scale_colour_manual(values = celltype_colours, name = "Celltype") +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  theme_classic()
ggsave(paste0(opts$output.base, "_joint_celltype.", opts$plot.type), plot = joint_celltype, device=opts$plot.type, width = 10, height = 5, units = "in")


joint_map <- ggplot(data=plot_df_joined, mapping = aes(x=X1, y=X2, colour=MAP)) +
  geom_point(size=io$dot_size, alpha=io$dot_alpha) +
  # ggrastr::geom_point_rast(size=opts$dot_size, alpha=opts$dot_alpha) +
  labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
  theme_classic()
ggsave(paste0(opts$output.base, "_joint_map.", opts$plot.type), plot = joint_map, device=opts$plot.type, width = 7, height = 5, units = "in")


to.plot <- meta_mapping[,.N,by=c("celltype.mapped","batch")]
celltypes_mapped <- barplot.pub(to.plot, x="celltype.mapped", colors=celltype_colours) +
  facet_wrap(~batch, nrow=1, scales="free_x")
ggsave(paste0(opts$output.base, "_celltypes_mapped.", opts$plot.type), plot = celltypes_mapped, device=opts$plot.type, width = 7, height = 5, units = "in")


to.plot <- meta_mapping
celltypes_score <- gghistogram(to.plot, x="celltype.score", y="..density..", facet="batch",
                               panel.labs = list(batch = unlist(lapply(sort(unique(to.plot$batch)), function(x) paste0(toString(x), ", n=", toString(nrow(to.plot[which(to.plot$batch==x),]))))))) +
  labs(x="Celltype mapping score (per query cell)", y="Density") +
  theme(legend.title = element_blank())
ggsave(paste0(opts$output.base, "_celltypes_score.", opts$plot.type), plot = celltypes_score, device=opts$plot.type, width = 5, height = 5, units = "in")


celltypes <- unique(meta_mapping$celltype.mapped)
for (ct in celltypes) {
    if (is.na(ct)) {
        to.plot <- meta_mapping[which(is.na(meta_mapping$celltype.mapped)),]
    } else {
        to.plot <- meta_mapping[which(meta_mapping$celltype.mapped==ct),]
    }

    p <- gghistogram(to.plot, x="celltype.score", y="..density..", facet="batch",
                     title=paste0(toString(ct),", n=",toString(nrow(to.plot))),
                     panel.labs = list(batch = unlist(lapply(sort(unique(to.plot$batch)), function(x) paste0(toString(x), ", n=", toString(nrow(to.plot[which(to.plot$batch==x),]))))))) +
      labs(x="Celltype mapping score (per query cell)", y="Density") +
      theme(legend.title = element_blank())
    ggsave(paste0(opts$output.base, "_celltypes_score_", gsub("/", "", toString(ct)),".", opts$plot.type), plot = p, device=opts$plot.type, width = 5, height = 5, units = "in")
}

to.plot <- meta_mapping
cellstages_score <- gghistogram(to.plot, x="cellstage.score", y="..density..", facet="batch",
            title=paste0(toString(ct),", n=",toString(nrow(to.plot))),
            panel.labs = list(batch = unlist(lapply(sort(unique(to.plot$batch)), function(x) paste0(toString(x), ", n=", toString(nrow(to.plot[which(to.plot$batch==x),]))))))) +
  labs(x="Cellstage mapping score (per query cell)", y="Density") +
  theme(legend.title = element_blank())
ggsave(paste0(opts$output.base, "_cellstages_score.", opts$plot.type), plot = cellstages_score, device=opts$plot.type, width = 5, height = 5, units = "in")
