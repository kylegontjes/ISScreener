# ISScreener
Runs insertion sequence detection workflow on Illumina WGS data

# Sample command
snakemake -s ISScreener.smk --use-conda --use-singularity -j 999 --cluster "sbatch -A {cluster.account} -p {cluster.partition} -N {cluster.nodes}  -t {cluster.walltime} -c {cluster.procs} --mem-per-cpu {cluster.pmem}" --conda-frontend conda --cluster-config config/cluster.json --configfile config/config.yaml --latency-wait 1000

# To run many isolates at the same time (and possibly close the computer, try running this command)
srun --account=esnitkin1 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=5GB --cpus-per-task=1 --time=12:00:00 --pty /bin/bash

# Curating the samples_list.csv file
## Given a path to a directory with files format

path="/nfs/esnitkin/Project_Penn_KPC/Sequence_data/fastq/Penn/SRA_submission/"

sample_id="sample_id"

sample_names=$(ls -1 $path | grep _R1 |  cut -d. -f1 | sed 's\_R1\\' | sed 's\_R2\\')

echo -e\n $sample_id $sample_names | tr ' ' '\n' > config/sample.tsv
