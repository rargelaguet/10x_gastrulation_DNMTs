source(here::here("settings.R"))
source(here::here("utils.R"))
source(here::here("pseudobulk/utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$repeats_expr <- file.path(io$basedir,"results/repeats/repeats_expr.txt.gz")
io$outdir <- file.path(io$basedir,"shiny/repeats")

opts$rename_celltypes <- c(
  "Erythroid3" = "Erythroid",
  "Erythroid2" = "Erythroid",
  "Erythroid1" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors"
)

##########################
## Load expression data ##
##########################

repeats.dt <- fread(io$repeats_expr) %>%
  .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,opts$rename_celltypes)] %>%
  .[,.(expr=mean(expr)),by=c("alias","celltype.mapped","class","repeat_class")] %>%
  setnames(c("sample","celltype","class","repeat_class","expr"))

# print
unique(repeats.dt$celltype)
unique(repeats.dt$repeat_class)

#######################################
## Calculate differential expression ##
#######################################

diff_repeats.dt <- repeats.dt %>% 
  .[,.(expr=mean(expr)),by=c("celltype","class","repeat_class")] %>% 
  dcast(celltype + repeat_class ~ class, value.var = "expr") %>% 
  .[, Dnmt1_KO := Dnmt1_KO - WT] %>% 
  .[, Dnmt3a_KO := Dnmt3a_KO - WT] %>% 
  .[, Dnmt3b_KO := Dnmt3b_KO - WT] %>% 
  melt(id.vars = c("celltype", "repeat_class"), measure.vars = c("Dnmt1_KO", "Dnmt3a_KO", "Dnmt3b_KO"), value.name = "logFC",variable.name = "class") %>%
  .[,logFC:=round(logFC,2)]

##########
## Save ##
##########

fwrite(repeats.dt, file.path(io$outdir,"repeats_expr.txt.gz"), na="NA", quote=F, sep="\t")
fwrite(diff_repeats.dt, file.path(io$outdir,"repeats_diff_expr.txt.gz"), na="NA", quote=F, sep="\t")
