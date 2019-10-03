#!/bin/bash
source activate v3
snakemake --keep-going -j 50 -s ./tiger.snek  --configfile config/config_rackham_all_samples.yaml  --cluster-config config/slurm.json  --cluster config/slurm_scheduler.py  --cluster-status config/slurm_status.py --rerun-incomplete --use-conda
