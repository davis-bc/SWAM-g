# EPA SWAM WGS
Snakemake runs on the node you launched it from (usually login/head or an interactive allocation).

Snakemake itself submits each rule as its own sbatch job, with its own SLURM header (taken from cluster.yaml + --cluster template).

The rules run on the cluster as independent jobs.

Downside: Snakemake stays alive on the head/login node while jobs run. Some HPCs discourage long-running processes on login nodes.

```bash
#!/bin/bash
#SBATCH --account=nrsaamr
#SBATCH --partition=ord
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2g
#SBATCH --job-name=snakemake_driver
#SBATCH --output=slurm_logs/snakemake.%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=davis.benjamin@epa.gov

module load snakemake/ # or source your env
snakemake --use-conda --jobs 50 \
          --config organism="Salmonella" \
          --cluster-config cluster.yaml \
          --cluster "sbatch --account={cluster.account} \
                            --partition={cluster.partition} \
                            --time={cluster.time} \
                            --mem={cluster.mem} \
                            --cpus-per-task={cluster.cpus-per-task} \
                            --mail-type={cluster.mail-type} \
                            --mail-user={cluster.mail-user} \
                            --output={cluster.output}"
```

sbatch run_snakemake.sh
