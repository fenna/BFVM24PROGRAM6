snakemake --snakefile main.smk --cores 4
snakemake --snakefile main.smk --dag | dot -Tpng > dag.png