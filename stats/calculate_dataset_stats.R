here::i_am("processing/stats/calculate_dataset_stats.R")

source(here::here("settings.R"))

#####################
## Define settings ##
#####################

# I/O
io$metadata <- file.path(io$basedir,"results_all/sex_assignment/sample_metadata_after_sex_assignment.txt.gz")
io$outdir <- file.path(io$basedir,"results_all/dataset_stats"); dir.create(io$outdir, showWarnings = F)

# Options
opts$classes <- c(
    "WT",
    "Dnmt3a_KO",
    "Dnmt3b_KO",
    "Dnmt1_KO",
    "Dnmt3ab_KO"
)

###############
## Load data ##
###############

metadata <- fread(io$metadata) %>% 
    .[,dataset:=ifelse(grepl("Grosswendt",sample),"CRISPR (Grosswendt2020)","KO")] %>%
    .[pass_rnaQC==TRUE & class%in%opts$classes] %>%
    .[,class:=factor(class, levels=opts$classes)]

#################################################
## Plot number of cells per class and data set ##
#################################################

to.plot <- metadata %>% .[,.N,by=c("class","dataset")]

p <- ggbarplot(to.plot, x="class", y="N", fill="dataset", position=position_dodge(width = 0.75)) +
    labs(x="", y="Number of cells") +
    scale_fill_brewer(palette="Dark2") +
    theme(
        legend.position = "top",
        legend.title = element_blank(),
        axis.text.y = element_text(colour="black",size=rel(0.8)),
        axis.text.x = element_text(colour="black",size=rel(0.80)),
    )

pdf(file.path(io$outdir,"ncells_per_class.pdf"), width=6, height=5.5)
print(p)
dev.off()

###################################################
## Plot number of embryos per class and data set ##
###################################################

to.plot <- metadata %>% .[,.(N=as.factor(length(unique(alias)))),by=c("class","dataset")]

p <- ggbarplot(to.plot, x="class", y="N", fill="dataset", position=position_dodge(width = 0.75)) +
    labs(x="", y="Number of embryos") +
    scale_fill_brewer(palette="Dark2") +
    theme(
        legend.position = "top",
        legend.title = element_blank(),
        axis.text.y = element_text(colour="black",size=rel(0.8)),
        axis.text.x = element_text(colour="black",size=rel(0.80)),
    )

pdf(file.path(io$outdir,"nembryos_per_class.pdf"), width=6, height=5.5)
print(p)
dev.off()

####################################################
## Plot number of embryos per class, split by sex ##
####################################################

to.plot <- metadata %>% .[,.(N=as.factor(length(unique(alias)))),by=c("class","sex")]

p <- ggbarplot(to.plot, x="class", y="N", fill="sex", position=position_dodge(width = 0.75)) +
    labs(x="", y="Number of embryos") +
    scale_fill_brewer(palette="Accent") +
    theme(
        legend.position = "top",
        legend.title = element_blank(),
        axis.text.y = element_text(colour="black",size=rel(0.8)),
        axis.text.x = element_text(colour="black",size=rel(0.80)),
    )

pdf(file.path(io$outdir,"nembryos_per_sex.pdf"), width=6, height=5.5)
print(p)
dev.off()
