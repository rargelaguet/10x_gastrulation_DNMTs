library(ggplot2)
linMap <- function(x, from, to) return( (x - min(x)) / max(x - min(x)) * (to - from) + from )

source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")
# source("/Users/ricard/10x_gastrulation_DNMTs/mapping/plot/plot_settings.R")

################
## Define I/O ##
################

io$outdir <- paste0(io$basedir,"/results/mapping/pdf")

####################
## Define options ##
####################

# Figure dimensions
# opts$width <- 3
# opts$height <- 4

opts$batches <- c(
  "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001", 
  "SIGAB6_E85_3_Dnmt3aWT_Dnmt3b_WT_L002", 
  "SIGAC6_E85_5_Dnmt3aKO_Dnmt3b_Het_L003", 
  "SIGAD6_E85_8_Dnmt3aHet_Dnmt3b_KO_L004",
  "15_E8_5_D3A_WT_D3B_WT_L007",
  "17_E8_5_D3A_KO_D3B_WT_L008",
  "2_E8_5_D3A_WT_D3B_KO_L003",
  "3_E8_5_D3A_HET_D3B_WT_L004",
  "7_E8_5_D3A_WT_D3B_KO_L005",
  "8_E8_5_D3A_KO_D3B_KO_L006",
  
  # this should behave as WT...
  "3_E8_5_D3A_HET_D3B_WT_L004"
  
  # "E125_DNMT3A_HET_A_L001", 
  # "E125_DNMT3A_HET_A_L003", 
  # "E125_DNMT3A_KO_B_L002", 
  # "E125_DNMT3A_KO_E_L004", 
  # "A_E12_5_D3a_Het_L001",
  # "B_E12_5_D3a_KO_L002"
)

opts$wt.batches <- c(
  "15_E8_5_D3A_WT_D3B_WT_L007",
  "SIGAB6_E85_3_Dnmt3aWT_Dnmt3b_WT_L002"
)
  
opts$to.merge <- c(
  "Erythroid3" = "Erythroid",
  "Erythroid2" = "Erythroid",
  "Erythroid1" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors",
  "Intermediate_mesoderm" = "Mixed_mesoderm",
  "Paraxial_mesoderm" = "Mixed_mesoderm",
  "Nascent_mesoderm" = "Mixed_mesoderm",
  "Pharyngeal_mesoderm" = "Mixed_mesoderm",
  "Visceral_endoderm" = "ExE_endoderm"
)

# Update sample metadata
sample_metadata <- sample_metadata %>%
  .[batch%in%opts$batches] %>%
  .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,opts$to.merge)]

# Calculate WT proportions
wt.dt <- sample_metadata[batch%in%opts$wt.batches] %>%
  .[,ncells:=.N, by="batch"] %>%
  # .[,.(proportion=.N/unique(ncells)),by=c("celltype.mapped","batch")] %>%
  .[,.(proportion=.N/unique(ncells), N=.N),by=c("celltype.mapped","batch")] %>%
  .[,.(proportion=mean(proportion), N=mean(N)),by="celltype.mapped"]
  # .[,c("class","batch"):="WT"]


#############################################
## Plot number of cells for each cell type ##
#############################################

ko.dt <- sample_metadata[!batch%in%opts$wt.batches] %>%
  .[,ncells:=.N, by="batch"] %>%
  .[,.(proportion=.N/unique(ncells), N=.N),by=c("celltype.mapped","batch","class")]
  

# Merge
dt <- merge(ko.dt, wt.dt, by="celltype.mapped", allow.cartesian=T, suffixes = c(".ko",".wt"))

# to.plot <- dt[,.(diff_proportion=log2(proportion.ko/proportion.wt), diff_N=abs(N.ko-N.wt)), by=c("batch","class","celltype.mapped")]# %>%
to.plot <- dt[,.(diff_proportion=log2(proportion.ko/proportion.wt), diff_N=N.ko-N.wt), by=c("batch","class","celltype.mapped")]# %>%
  # setnames(c("batch","class","celltype.mapped","diff"))

to.plot[,diff_N_norm:=linMap(abs(diff_N), from=0.15, to=2.5)]

# sequence_length = length(unique(to.plot.test$celltype.mapped))
# first_sequence = c(1:(sequence_length%/%2)) 
# second_sequence = c((sequence_length%/%2+1):sequence_length) 
# first_angles = c(90 - 180/length(first_sequence) * first_sequence)
# second_angles = c(-90 - 180/length(second_sequence) * second_sequence)

to.plot.wt_line <- data.table(
  celltype.mapped = unique(dt$celltype.mapped),
  diff_proportion = log2(1),
  diff_N = 0 #+  0.01
)


ggplot(to.plot, aes(x=factor(celltype.mapped), y=diff_proportion, group=1)) +
  geom_point(aes(color = celltype.mapped, size = diff_N_norm), stat = 'identity') + 
  geom_polygon(color="black", fill=NA, alpha=0.5, linetype="dashed", data=to.plot.wt_line) +
  facet_wrap(~batch) +
  scale_color_manual(values=opts$celltype.colors, drop=F) +
  coord_polar() +
  guides(colour = guide_legend(override.aes = list(size=2), ncol=1)) +
  scale_size(guide = 'none') +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.text = element_text(size=rel(0.75)),
    legend.title = element_blank(),
    axis.title=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.line=element_blank(),
    axis.text.x = element_blank(),
    # axis.text.x = element_text(angle= -76 - 360 / length(unique(to.plot.test$celltype.mapped)) * seq_along(to.plot.test$celltype.mapped))
  )

# pdf(paste0(io$outdir,"/pdf/mapping_stats_day2.pdf"), width=opts$width, height=opts$height)
# print(p)
# dev.off()






########################################
## Plot difference in number of cells ##
########################################

# ggplot(to.plot, aes(x=factor(celltype.mapped), y=log2(diff_N+0.01), group=1)) +
ggplot(to.plot, aes(x=factor(celltype.mapped), y=diff_N, group=1)) +
  geom_point(aes(color = celltype.mapped, size = diff_N_norm), stat = 'identity') + 
  geom_polygon(color="black", fill=NA, alpha=0.5, linetype="dashed", data=to.plot.wt_line) +
  facet_wrap(~batch) +
  scale_color_manual(values=opts$celltype.colors, drop=F) +
  coord_polar() +
  guides(colour = guide_legend(override.aes = list(size=2), ncol=1)) +
  scale_size(guide = 'none') +
  theme_bw() +
  theme(
    legend.position = "none",
    legend.text = element_text(size=rel(0.75)),
    legend.title = element_blank(),
    axis.title=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.line=element_blank(),
    axis.text.x = element_blank(),
    # axis.text.x = element_text(angle= -76 - 360 / length(unique(to.plot.test$celltype.mapped)) * seq_along(to.plot.test$celltype.mapped))
  )
