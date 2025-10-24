#!/bin/bash
#SBATCH --account=def-calcium
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4096M
#SBATCH --time=0-2:00
#SBATCH --mail-user=qihuang.zhang@mcgill.ca
#SBATCH --mail-type=ALL

cd /lustre03/project/6075067/calcium/2022/GSIMEX
module load r/4.2.2

R CMD BATCH '--args range=c(301, 400)' code/Simulations/LM/GSIMEX-LM.R code/Simulations/LM/GSIMEX-LM.out 