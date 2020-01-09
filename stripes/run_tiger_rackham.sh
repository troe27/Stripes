#!/bin/bash
source activate v3
snakemake --keep-going -j 150 -s ./tiger.snek --configfile config/config.yaml  --cluster-config config/slurm.json --unlock  --cluster config/slurm_scheduler.py  --cluster-status config/slurm_status.py --rerun-incomplete --use-conda
