# Get the directory of the current R script
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

# Set the working directory to the directory of the script
setwd(script_dir)

for (i in 1:10) {
  lower_range <- 1 + 50 * (i - 1)
  upper_range <- 50 * i
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

R CMD BATCH '--args range=c(%d, %d)' code/Simulations/GSIMEX-ZIP.R code/Simulations/GSIMEX-ZIP.out ",
    lower_range, upper_range
  )
  
  script_filename <- sprintf("GSIMEX-ZIP%d.sh", i)
  
  # Write the script to a file
  cat(script, file = script_filename)
  
  # Print a message indicating the script generation
  cat(sprintf("Generated script: %s\n", script_filename))
}