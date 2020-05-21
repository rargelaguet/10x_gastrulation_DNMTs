# mapping_TetTKO
Git repo with pipeline for mapping Tet TKO chimeras, given an already processed Seurat.

NOTE: Some files are named DNMT3A when in fact they're DNMT3B. I have renamed this in the metadata file when parsing the modifications and targets, but kept the batch name, as for better or worse this is the file name.

This pipeline was built with Snakemake (https://snakemake.readthedocs.io/en/stable/) and so it can be run from the EBI cluster with the following command: - bsub -M 50 snakemake -j 99 --cluster "bsub {params.other} -M {params.memory} -R rusage[mem={params.rusage}] -n 1 -o {params.outfile} -e {params.errfile} -q research-rh74".

## Editing
To rerun the pipeline, the user need only edit the config.yaml file with the appropriate directories and file paths.
The pipeline was built for an lsf cluster, and so if your cluster does not use lsf, you may need to change some params in the snakemake rules appropriately.
Each step has a settings.R file to go with it, which gives the option of hard-coding inputs/outputs, and which also hard-codes samples and sample-specific information. If this pipeline is run on new data, these will have to be ammended.

## File Structure
There are 3 steps:
1. pre-processing (this takes the Seurat object and does some further QC)
2. mapping (this does the actual mapping). Within the mapping directory, there is a separate directory for each method (harmony,mnn,seurat,liger). Harmony, MNN, and LIGER map the atlas new each time, because they are very fast. Seurat takes a lot longer to map the atlas, so does it once, and this mapped atlas is then input for each experiment individually.
3. evaluating. Here the output plots are made.

The mapping and evaluating steps run twice: Once on the whole dataset, and a second time it subsets cells based on their previous mapping and runs on each subset individually. The subsets are: "Haematoendothelial", "Endoderm", "EpiblastPSNeuro", and "Mesoderm" (set in /mapping/settings.R).
    
