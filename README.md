# ISScreener
Runs insertion sequence detection workflow on Illumina WGS data

# Sample command
snakemake -s ISScreener.smk -p --use-conda --use-singularity -j 999 --cluster "sbatch -A {cluster.account} -p {cluster.partition} -N {cluster.nodes}  -t {cluster.walltime} -c {cluster.procs} --mem-per-cpu {cluster.pmem}" --conda-frontend conda --cluster-config config/cluster.json --configfile config/config.yaml --latency-wait 1000
