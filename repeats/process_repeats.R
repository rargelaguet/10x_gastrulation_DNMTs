# here::i_am("processing/1_create_seurat_rna.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

# suppressPackageStartupMessages(library(Seurat))

#####################
## Define settings ##
#####################

io$inputdir <- file.path(io$basedir,"processed/repeats")
io$outdir <- file.path(io$basedir,"results/repeats"); dir.create(io$outdir, showWarnings = F)

opts$sample2class <- c(
	Dnmt1_KO_1 = "Dnmt1_KO",
	Dnmt1_KO_2 = "Dnmt1_KO",
	Dnmt1_KO_3 = "Dnmt1_KO",
	Dnmt3a_KO_13 = "Dnmt3a_KO",
	Dnmt3a_KO_14 = "Dnmt3a_KO",
	Dnmt3a_KO_1 = "Dnmt3a_KO",
	Dnmt3a_KO_2 = "Dnmt3a_KO",
	Dnmt3b_KO_11 = "Dnmt3b_KO",
	Dnmt3b_KO_12 = "Dnmt3b_KO",
	Dnmt3b_KO_1 = "Dnmt3b_KO",
	Dnmt3b_KO_2 = "Dnmt3b_KO",
	WT_1 = "WT",
	WT_2 = "WT",
	WT_3 = "WT",
	WT_4 = "WT",
	WT_5 = "WT",
	WT_6 = "WT",
	WT_7 = "WT"
)

repeats_names <- c("LINE_L1", "LINE_L2", "LTR_ERV1", "LTR_ERVK", "LTR_ERVL", "LTR_MaLR", "major_satellite", "minor_satellite", "rRNA", "SINE_Alu_B1", "SINE_B2", "SINE_B4", "telomere", "IAP")

###################
## Load metadata ##
###################

cell_metadata.dt <- fread(io$metadata) %>% 
  .[pass_rnaQC==TRUE & alias%in%names(opts$sample2class)]
table(cell_metadata.dt$alias)

cell_metadata.dt[,barcode:=gsub("-1","",barcode)]

#####################################
## Load and merge repeat data sets ##
#####################################

# i <- "WT_1"
repeats_wide.dt <- names(opts$sample2class) %>% map(function(i) {
	tmp <- fread(file.path(io$inputdir,sprintf("%s.txt.gz",i)), sep="\t") %>% 
	  .[,barcode:=gsub("CB:","",barcode)] %>%
	  .[barcode%in%cell_metadata.dt[alias==i,barcode]] %>%
	  .[,alias:=i] %>% .[,foo:=NULL]
}) %>% rbindlist

cell_metadata.keys <- c("barcode","alias","celltype.mapped","dataset","class","celltype_class")

total_reads.dt <- repeats_wide.dt[,c("barcode","total_reads","alias")] %>% .[,.(total_reads=sum(total_reads)),by="alias"]

repeats_long.dt <- merge(repeats_wide.dt,cell_metadata.dt[,cell_metadata.keys,with=F],by=c("barcode","alias")) %>% 
  melt(id.vars=c("barcode","alias","celltype.mapped","dataset","class","celltype_class"), measure.vars=repeats_names, variable.name="repeat_class", value.name="nreads") %>%
  .[,.(nreads=sum(nreads)),by=c(cell_metadata.keys[2:length(cell_metadata.keys)],"repeat_class")]


###############
## Normalise ##
###############

repeats_long.dt <- repeats_long.dt %>% 
  merge(total_reads.dt,by=c("alias")) %>% 
  .[,expr:=log2(1+(1e6*(nreads/total_reads))) %>% round(2)]

######################
## Compare WT vs KO ##
######################

tmp <- cell_metadata.dt[,.N,by=c("class","celltype.mapped")] %>% .[N>=100]

celltypes.to.plot <- tmp[,.N,by="celltype.mapped"] %>% .[N>=3,celltype.mapped] %>% unique

# i <- "LINE_L1"
for (i in repeats_names) {
  
  to.plot <- repeats_long.dt[repeat_class==i] %>% 
    .[celltype.mapped%in%celltypes.to.plot] %>%
    merge(tmp,by=c("class","celltype.mapped"))
  
  p <- ggplot(to.plot, aes(x=class, y=expr, fill=celltype.mapped)) +
    geom_boxplot(aes(fill = class), alpha=0.75) +
    facet_wrap(~celltype.mapped, scales="free_y") +
    scale_fill_manual(values=opts$classes.colors[unique(to.plot$class)]) +
    stat_summary(fun.data = give.n, geom = "text", size=3) +
    # stat_compare_means(aes(label = paste0("p = ", ..p.format..)), method="t.test") +
    # stat_compare_means(label="p.signif", label.y = 70, hide.ns=T)
    theme_classic() +
    labs(y="", x="") +
    theme(
      legend.position = "right",
      axis.text.y = element_text(color="black", size=rel(0.8)),
      axis.text.x = element_blank()
    )
  
  # pdf(file.path(io$outdir,sprintf("%s_boxplots.pdf",i)), width=6, height=7)
  pdf(file.path(io$outdir,sprintf("%s_repeats_boxplots.pdf",i)), width=12, height=7)
  print(p)
  dev.off()
}




##############################
## Compare embryonic vs ExE ##
##############################

i <- "LINE_L1"
to.plot <- repeats_long.dt[repeat_class==i] %>% 
  .[celltype.mapped%in%celltypes.to.plot] %>%
  merge(tmp,by=c("class","celltype.mapped"))

p <- ggplot(to.plot, aes(x=celltype.mapped, y=expr)) +
  geom_boxplot(aes(fill=celltype.mapped), alpha=0.75) +
  scale_fill_manual(values=opts$celltype.colors[celltypes.to.plot]) +
  facet_wrap(~class, scales="fixed", nrow=1) +
  stat_summary(fun.data = give.n, geom = "text", size=3) +
  theme_bw() +
  labs(y="", x="RNA expression") +
  theme(
    legend.position = "none",
    axis.text.y = element_text(color="black", size=rel(0.8)),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# pdf(file.path(io$outdir,sprintf("%s_boxplots.pdf",i)), width=6, height=7)
pdf(file.path(io$outdir,sprintf("%s_repeats_boxplots_v2.pdf",i)), width=12, height=7)
print(p)
dev.off()