#!/bin/bash
#SBATCH -A m342               # project account
#SBATCH -C cpu                # request CPU‐only nodes
#SBATCH -q regular            # shared queue
#SBATCH -t 10:59:00           # walltime (HH:MM:SS)
#SBATCH -n 1                  # number of tasks
#SBATCH -c 256                # number of CPU cores per task
#SBATCH -J grassfam           # job name
#SBATCH -o logs/%x-%j.out     # STDOUT log (%x = job name, %j = job ID)
#SBATCH -e logs/%x-%j.err     # STDERR log

srun bash main.sh 250 &> log.log