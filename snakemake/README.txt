snakemake -j 1
bsub -M 50 snakemake -j 1 --cluster "bsub -M {params.memory} -n 1 -q research-rh74"
bsub -M 50 snakemake -j 99 --cluster "bsub {params.other} -M {params.memory} -R rusage[mem={params.rusage}] -n 1 -o {params.outfile} -e {params.errfile} -q research-rh74"

TO-DO:
- soupX
- velocity



