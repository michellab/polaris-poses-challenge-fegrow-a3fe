#!/bin/bash
#SBATCH -o somd-array-gpu-%A.%a.out
#SBATCH -p RTX4060,RTX3080
#SBATCH -n 1
#SBATCH --time 24:00:00
#SBATCH --gres=gpu:1

lam=$1
echo "lambda is: " $lam

srun somd-freenrg -C somd.cfg -l $lam -p CUDA

