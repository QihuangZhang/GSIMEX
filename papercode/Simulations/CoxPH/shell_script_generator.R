# Get the directory of the current R script
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

# Set the working directory to the directory of the script
setwd(script_dir)

node_n <- 5
job_per_node <- 100

for (i in 1:node_n) {
  lower_range <- 1 + job_per_node * (i - 1)
  upper_range <- job_per_node * i
  script <- sprintf(
    "#!/bin/bash
#SBATCH --account=def-calcium
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4096M
#SBATCH --time=0-2:00
#SBATCH --mail-user=qihuang.zhang@mcgill.ca
#SBATCH --mail-type=ALL

cd /lustre03/project/6075067/calcium/2022/GSIMEX
module load r/4.2.2

R CMD BATCH '--args range=c(%d, %d)' code/Simulations/coxPH/GSIMEX-coxPH.R code/Simulations/coxPH/GSIMEX-coxPH.out ",
    lower_range, upper_range
  )
  
  script_filename <- sprintf("GSIMEX-coxPH%d.sh", i)
  
  # Write the script to a file
  cat(script, file = script_filename)
  
  # Print a message indicating the script generation
  cat(sprintf("Generated script: %s\n", script_filename))
}