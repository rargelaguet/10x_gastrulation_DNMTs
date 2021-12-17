#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(MouseGastrulationData))
suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(batchelor))
#suppressPackageStartupMessages(library(future))
#plan("multiprocess", workers=2) # parallelise to use all available cores
options(future.globals.maxSize= 115000*1024^2)


###################
## Parse options ##
###################
 
option_list = list(
    make_option(c("-A", "--atlas.RDS"), type="character", default=NULL, 
                help="path to the mouse gastrulation 10x atlas SCE RDS", metavar="character"),
    make_option(c("-a", "--atlas.metadata"), type="character", default=NULL, 
                help="path to the mouse gastrulation 10x atlas metdata txt", metavar="character"),
    make_option(c("-O", "--output.RDS"), type="character", default=NULL, 
                help="path to output RDS file", metavar="character"),
    make_option(c("-o", "--output.metadata"), type="character", default=NULL, 
                help="path to output sample metadata file", metavar="character"),
    make_option(c("-l", "--lineage.list"), type="character", default=NULL, 
                help="which lineages to use", metavar="character"),
    make_option(c("-L", "--lineage.name"), type="character", default=NULL, 
                help="name of the lineage subset to use", metavar="character"),
    make_option(c("-s", "--atlas_stages"), type="character", default=NULL, 
                help="stages of the atlas to map to", metavar="character"),
    make_option(c("-S", "--settings"), type="character", default=NULL, 
                help="path to settings file", metavar="character"),
    make_option(c("-b", "--query_batches"), type="character", default=NULL, 
                help="names of batches in query", metavar="character"),
    make_option(c("-t", "--testing"), action="store_true", default=NULL, 
                help="if given, data will be subset, so code can be tested")
); 
 
opt_parser = OptionParser(option_list=option_list);
opts = parse_args(opt_parser);

if (is.null(opts$settings)) {
    print_help(opt_parser)
    stop("A settings file must be supplied.n", call.=FALSE)
}

source(opts$settings)

if ((!is.null(io$testing)) && (io$testing != TRUE)){
    stop("In the settings file, testing must be NULL or TRUE.n", call.=FALSE)
} else if ((!is.null(io$atlas.metadata)) && (!is.string(io$atlas.metadata))){
    stop("In the settings file, atlas.metadata must be NULL or a path.n", call.=FALSE)
} else if ((!is.null(io$atlas.RDS)) && (!is.string(io$atlas.RDS))){
    stop("In the settings file, atlas.RDS must be NULL or a path.n", call.=FALSE)
} else if ((!is.null(io$atlas_stages)) && (!is.string(io$atlas_stages))){
    stop("In the settings file, atlas_stages must be NULL or a path.n", call.=FALSE)
} else if ((!is.null(opts$lineage.name)) && (!opts$lineage.name %in% names(io$subset_celltypes))){
    stop("The lineage name must be the name of one of the settings file classes.n", call.=FALSE)
}

if (isTRUE(io$testing)) opts$testing = TRUE
if (is.string(io$atlas.RDS)) opts$atlas.RDS <- io$atlas.RDS
if (is.string(io$atlas.metadata)) opts$atlas.metadata <- io$atlas.metadata
if (is.string(io$atlas_stages)) opts$atlas_stages <- io$atlas_stages
if (isTRUE(opts$testing)) message("Doing mapping in test mode (subsetting cells)...")

if (is.null(opts$atlas.RDS)) {
    print_help(opt_parser)
    stop("A path to an RDS input file must be supplied.n", call.=FALSE)
} else if (is.null(opts$atlas.metadata)) {
    print_help(opt_parser)
    stop("A path to a metadata input file must be supplied.n", call.=FALSE)
} else if (is.null(opts$output.RDS) && !isTRUE(opts$testing)) {
    print_help(opt_parser)
    stop("A path to an RDS output file must be supplied.n", call.=FALSE)
} else if (is.null(opts$output.metadata) && !isTRUE(opts$testing)) {
    print_help(opt_parser)
    stop("A path to a metadata output file must be supplied.n", call.=FALSE)
}

if (is.string(opts$lineage.list)) opts$lineage.list <- str_trim(strsplit(opts$lineage.list, ",")[[1]])
if (is.string(opts$lineage.name)) opts$lineage.list <- io$subset_celltypes[opts$lineage.name][[1]]

################
## Load atlas ##
################

message("Loading atlas...")

if ((is.null(opts$atlas.RDS)) || (is.null(opts$atlas.metadata))){
    
    message("Downloading atlas from bioconda package...")
    
    # Load SingleCellExperiment
    suppressMessages(sce_atlas <- EmbryoAtlasData(type="processed"))
    
    sce_atlas <- sce_atlas[,sce_atlas$doublet==F & sce_atlas$stripped==F]
    
    # Load metadata
    meta_atlas <- colData(sce_atlas) %>% as.data.table
    
} else {
    
    # Load metadata
    if (file_ext(opts$atlas.metadata) == "rds") {
        meta_atlas <- readRDS(paste0(opts$atlas.metadata)) %>% as.data.table
    } else {
        meta_atlas <- fread(opts$atlas.metadata, sep="\t", na="NA", quote=F) %>% as.data.table
    }
    
    meta_atlas <- meta_atlas[meta_atlas$doublet==FALSE,]
    meta_atlas <- meta_atlas[meta_atlas$stripped==FALSE,]
    
    # Load SingleCellExperiment
    sce_atlas  <- readRDS(paste0(opts$atlas.RDS))[,meta_atlas$cell]
}

# Filter stages
if (!is.null(opts$atlas_stages)) {
  message(sprintf("Subsetting atlas stages to %s",paste(opts$atlas_stages, collapse=", ")))
  meta_atlas <- meta_atlas[stage%in%opts$atlas_stages]
}

# Filter lineages
if (!is.null(opts$lineage.list)) {
  message(sprintf("Subsetting atlas to %s",paste(opts$lineage.list, collapse=", ")))
  meta_atlas <- meta_atlas[meta_atlas$celltype %in% opts$lineage.list]
}
atlas_tpoints <- lapply(unique(meta_atlas$stage), toString)
atlas_samples <- lapply(unique(meta_atlas$sample), toString)
atlas_tpsamps <- unique(meta_atlas[,c("stage","sample")])

# Filter if test mode
if (isTRUE(opts$testing)) {
    
    n <- round(10000/length(atlas_samples))
    sub = c()
    for (s in atlas_samples){
        sub <- append(sub, meta_atlas[meta_atlas$sample==s,]$cell[1:n])
    }
    
    meta_atlas <- meta_atlas[meta_atlas$cell %in% sub,]
}


sce_atlas <- sce_atlas[,meta_atlas$cell]
stopifnot(all(meta_atlas$cell == colnames(sce_atlas)))

# Log Normalise
message("Normalising data...")
sce_atlas <- multiBatchNorm(sce_atlas, batch=meta_atlas$sample, preserve.single=TRUE)


###############
## Find HVGs ##
###############

message("Finding HVGs...")
decomp <- modelGeneVar(sce_atlas, block=meta_atlas$sample)
decomp <- decomp[decomp$mean > io$min.mean,]
decomp$FDR <- p.adjust(decomp$p.value, method = "fdr")
hvgs    <- rownames(decomp[order(decomp$p.value),][1:2000,])

# Create Seurat object
message("Creating Seurat object...")
seurat_atlas <- as.Seurat(x=sce_atlas, project='ATLAS')
rm(sce_atlas)
seurat_atlas@meta.data$map <- "ATLAS"

seurat_atlas@reductions$umap <- NULL
seurat_atlas@reductions$pca.corrected <- NULL
VariableFeatures(seurat_atlas) <- hvgs


#####################
## Split  by batch ##
#####################

message("Splitting by batch...")
seurat_atlas.list <- SplitObject(seurat_atlas, split.by = "sample")


#####################
## Integrate Atlas ##
#####################

message("Integrating datasets...")
seurat_atlas.sampintegrated <- list()
oldtp <- NULL
for (tp in atlas_tpoints) {
    sublist <- seurat_atlas.list[sapply(atlas_tpsamps[which(atlas_tpsamps$stage %in% c(oldtp, tp))]$sample,toString)]
    if (min(unlist(lapply(sublist, function(x) dim(x)[2]))) < 30) {
        tmp <- subset(seurat_atlas, subset = stage==tp)
        if (dim(tmp)[2] < 30) {
            oldtp <- c(oldtp, tp)
        } else {
            seurat_atlas.sampintegrated[tp] <- tmp
            oldtp <- NULL
        }
    } else if (length(sublist) == 1){
        seurat_atlas.sampintegrated[tp] <- sublist[[1]]
        oldtp <- NULL
    } else {
        integration.anchors <- FindIntegrationAnchors(object.list=sublist, k.filter=min(200, unlist(lapply(sublist, function(x) dim(x)[2]))), verbose=FALSE)
        seurat_atlas.sampintegrated[tp] <- IntegrateData(anchorset = integration.anchors, verbose=FALSE)
        oldtp <- NULL
    }
    
}
rm(sublist); rm(seurat_atlas.list)

integration.anchors <- FindIntegrationAnchors(object.list=seurat_atlas.sampintegrated, verbose=FALSE, k.filter=min(200, unlist(lapply(seurat_atlas.sampintegrated, function(x) dim(x)[2]))))
seurat_atlas.integrated <- IntegrateData(anchorset = integration.anchors, verbose=FALSE, k.weight=min(100, unlist(lapply(seurat_atlas.sampintegrated, function(x) dim(x)[2]))))
rm(seurat_atlas.sampintegrated); rm(integration.anchors); rm(seurat_atlas)

#############################
## Process Integrated Data ##
#############################

message("Scaling Data...")
seurat_atlas.integrated <- ScaleData(object = seurat_atlas.integrated, verbose=FALSE)
message("Dimensionality reduction...")
seurat_atlas.integrated <- RunPCA(seurat_atlas.integrated, npcs = io$npcs, verbose=FALSE)
seurat_atlas.integrated <- RunUMAP(seurat_atlas.integrated, umap.method = "umap-learn", reduction = "pca", dims = 1:io$npcs, verbose=FALSE)


##########
## Save ##
##########

message("Saving output...")
if (!isTRUE(opts$testing)) {
    
    # save mapping in .rds format
    saveRDS(seurat_atlas.integrated, paste0(opts$output.RDS))
    
    # save mapping in .txt format
    to.save <- seurat_atlas.integrated@meta.data %>% as.data.table
    fwrite(to.save, paste0(opts$output.metadata), sep="\t", na="NA", quote=F)
    
}




