# EPA SWAM WGS
Example code to run 

```bash
snakemake --use-conda \
          --jobs 50 \
          --cluster-config cluster.yaml \
          --config organism="Salmonella" \
          --cluster "sbatch --account={cluster.account} \
                            --partition={cluster.partition} \
                            --time={cluster.time} \
                            --mem={cluster.mem} \
                            --cpus-per-task={cluster.cpus-per-task} \
                            --mail-type={cluster.mail-type} \
                            --mail-user={cluster.mail-user} \
                            --output={cluster.output}"
```
