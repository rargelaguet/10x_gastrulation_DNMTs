linMap <- function(x, from, to) return( (x - min(x)) / max(x - min(x)) * (to - from) + from )

source("/Users/ricard/10x_gastrulation_DNMTs/settings.R")

################
## Define I/O ##
################

io$outdir <- paste0(io$basedir,"/results/mapping/pdf")

####################
## Define options ##
####################

opts$classes <- c(
  # "E12.5_Dnmt3aWT_Dnmt3bHET",
  # "E12.5_Dnmt3aWT_Dnmt3bKO",
  # "E12.5_Dnmt3aHET_Dnmt3bWT",
  # "E12.5_Dnmt3aKO_Dnmt3bWT",
  "E8.5_Dnmt3aKO_Dnmt3bWT",
  "E8.5_WT",
  "E8.5_Dnmt3aHET_Dnmt3bKO",
  "E8.5_Dnmt3aHET_Dnmt3bWT",
  "E8.5_Dnmt3aKO_Dnmt3bHET",
  "E8.5_Dnmt3aKO_Dnmt3bKO",
  "E8.5_Dnmt3aWT_Dnmt3bKO",
  "E8.5_Dnmt1KO"
  )

opts$wt.classes <- c("E8.5_WT")

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

############################
## Update sample metadata ##
############################

sample_metadata <- fread(io$metadata) %>%
  .[pass_QC==TRUE] %>%
  .[class%in%opts$classes & !is.na(celltype.mapped)] %>%
  # .[!celltype.mapped%in%c("ExE_ectoderm","ExE_endoderm","Parietal_endoderm","Visceral_endoderm")] %>%
  .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,opts$to.merge)]

table(sample_metadata$batch)

#########################
# Perform calculations ##
#########################



#########################################################
## Plot differences in cell type proportions per batch ##
#########################################################

# Calculate background proportions
# background.dt <- sample_metadata %>%
#   .[class%in%opts$wt.classes] %>%
#   .[,ncells:=.N, by="batch"] %>%
#   .[,.(proportion=.N/unique(ncells), N=.N),by=c("celltype.mapped","class")]# %>%
# # .[,.(proportion=mean(proportion), N=mean(N)),by="celltype.mapped"]
background.dt <- sample_metadata %>%
  .[class%in%opts$wt.classes] %>%
  .[,ncells:=.N] %>%
  .[,.(proportion=.N/unique(ncells), N=.N),by=c("celltype.mapped")]# %>%

# Calculate proportions
foreground.dt <- sample_metadata %>%
  # .[!class%in%opts$wt.classes] %>%
  .[,ncells:=.N, by="batch"] %>%
  .[,.(proportion=.N/unique(ncells), N=.N),by=c("celltype.mapped","batch","class")]

# Merge
dt <- merge(foreground.dt, background.dt, by=c("celltype.mapped"), allow.cartesian=T, suffixes = c(".ko",".wt"))

# Plot
to.plot <- dt %>%
  .[N.ko+N.wt>25] %>% # only consider cell types with enough observations
  # .[,.(diff_proportion=log2(proportion.ko/proportion.wt), diff_N=N.ko-N.wt), by=c("batch","celltype.mapped","class")] %>% 
  .[,.(diff_proportion=log2(proportion.ko/proportion.wt), diff_N=N.ko-N.wt), by=c("celltype.mapped","batch","class")] %>% 
  .[,diff_N_norm:=linMap(abs(diff_N), from=0.15, to=1.5)]

to.plot.wt_line <- data.table(
  celltype.mapped = unique(dt$celltype.mapped),
  diff_proportion = log2(1),
  diff_N = 0 #+  0.01
)

p <- ggplot(to.plot, aes(x=factor(celltype.mapped), y=diff_proportion, group=1)) +
  # geom_point(aes(color = celltype.mapped, size = diff_N_norm), stat = 'identity') + 
  geom_point(aes(color = celltype.mapped), size=2.5, stat = 'identity') + 
  geom_polygon(color="black", fill=NA, alpha=0.5, linetype="dashed", data=to.plot.wt_line) +
  # facet_wrap(~batch) +
  facet_wrap(~batch) +
  scale_color_manual(values=opts$celltype.colors, drop=F) +
  coord_polar() +
  guides(colour = guide_legend(override.aes = list(size=2), ncol=1)) +
  scale_size(guide = 'none') +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text =  element_text(size=rel(0.75)),
    legend.text = element_text(size=rel(0.75)),
    legend.title = element_blank(),
    axis.title=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.line=element_blank(),
    axis.text.x = element_blank()
    # axis.text.x = element_text(angle= -76 - 360 / length(unique(to.plot.test$celltype.mapped)) * seq_along(to.plot.test$celltype.mapped))
  )
  
pdf(sprintf("%s/polar_plots_by_baatch.pdf",io$outdir))
print(p)
dev.off()


#########################################################
## Plot differences in cell type proportions per class ##
#########################################################

# Calculate background proportions
# background.dt <- sample_metadata %>%
#   .[class%in%opts$wt.classes] %>%
#   .[,ncells:=.N, by="batch"] %>%
#   .[,.(proportion=.N/unique(ncells), N=.N),by=c("celltype.mapped","class")]# %>%
# # .[,.(proportion=mean(proportion), N=mean(N)),by="celltype.mapped"]
background.dt <- sample_metadata %>%
  .[class%in%opts$wt.classes] %>%
  .[,ncells:=.N] %>%
  .[,.(proportion=.N/unique(ncells), N=.N),by=c("celltype.mapped")]# %>%

# Calculate proportions
foreground.dt <- sample_metadata %>%
  # .[!class%in%opts$wt.classes] %>%
  .[,ncells:=.N, by="class"] %>%
  .[,.(proportion=.N/unique(ncells), N=.N),by=c("celltype.mapped","class")]

# Merge
dt <- merge(foreground.dt, background.dt, by=c("celltype.mapped"), allow.cartesian=T, suffixes = c(".ko",".wt"))

# Plot
to.plot <- dt %>%
  .[N.ko+N.wt>25] %>% # only consider cell types with enough observations
  .[,.(diff_proportion=log2(proportion.ko/proportion.wt), diff_N=N.ko-N.wt), by=c("celltype.mapped","class")] %>% 
  .[,diff_N_norm:=linMap(abs(diff_N), from=0.15, to=1.5)]

to.plot.wt_line <- data.table(
  celltype.mapped = unique(dt$celltype.mapped),
  diff_proportion = log2(1),
  diff_N = 0 #+  0.01
)

p <- ggplot(to.plot, aes(x=factor(celltype.mapped), y=diff_proportion, group=1)) +
  # geom_point(aes(color = celltype.mapped, size = diff_N_norm), stat = 'identity') + 
  geom_point(aes(color = celltype.mapped), size=2.5, stat = 'identity') + 
  geom_polygon(color="black", fill=NA, alpha=0.5, linetype="dashed", data=to.plot.wt_line) +
  # facet_wrap(~batch) +
  facet_wrap(~class) +
  scale_color_manual(values=opts$celltype.colors, drop=F) +
  coord_polar() +
  guides(colour = guide_legend(override.aes = list(size=2), ncol=1)) +
  scale_size(guide = 'none') +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text =  element_text(size=rel(0.9)),
    legend.text = element_text(size=rel(0.75)),
    legend.title = element_blank(),
    axis.title=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.line=element_blank(),
    axis.text.x = element_blank()
    # axis.text.x = element_text(angle= -76 - 360 / length(unique(to.plot.test$celltype.mapped)) * seq_along(to.plot.test$celltype.mapped))
  )

pdf(sprintf("%s/polar_plots_by_class.pdf",io$outdir))
print(p)
dev.off()