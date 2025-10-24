# GSIMEX
Generalized SIMEX with polynomial extrapolation and best-subset selection.

This repository contains code to reproduce the simulations and the real-data application in our manuscript on **GSIMEX** (generalized SIMEX). The implementation is in R with a few helper shell scripts.

---

## Contents
```
GSIMEX/
├─ papercode/                 # all analysis code lives here
│  ├─ simulations/            # Simulation 1–3 scripts
│  ├─ application/            # real-data (snRNA-seq + CeLEry) scripts
│  ├─ utility.R               # helper functions (plotting, IO, GSIMEX helpers)
├─ output/                    # figures and intermediate results (created at runtime)
├─ data/                      # (you create this) input data or symlinks
└─ README.md
```
---

## Requirements

- **R** ≥ 4.2 (tested on 4.3/4.4)
- Recommended R packages (installed automatically below):
  - `tidyverse`, `data.table`, `ggplot2`, `patchwork`, `viridisLite`
  - `glmnet`, `MASS`, `survival` (for Cox PH sims)
  - `foreach`, `doParallel`, `optparse`
  - For the real-data app: `Seurat` (or `SeuratObject`), `Matrix`, and dependencies used by **CeLEry**
- **CeLEry** package (for layer inference): see its repository for installation instructions.
- Optional: `renv` or `pak` for reproducible package installs.

### Quick setup

```r
# inside R
pkgs <- c("tidyverse","data.table","ggplot2","patchwork","viridisLite",
          "glmnet","MASS","survival","foreach","doParallel","optparse",
          "Seurat","SeuratObject","Matrix")
install.packages(setdiff(pkgs, rownames(installed.packages())))
```

---

## Getting started

```bash
git clone https://github.com/QihuangZhang/GSIMEX.git
cd GSIMEX
mkdir -p output data
```

> If you use `renv`:
> ```r
> install.packages("renv"); renv::init(); renv::snapshot()
> ```

---

## Reproducing simulations

All scripts below write figures to `output/Simulation/` and save intermediate `.rds` files so you can re-run plotting quickly.

> **Note:** Script names below reflect the typical layout; if your repo uses different filenames, just swap them in the commands.

---

## Citation
If you use this code, please cite the manuscript and this repository.

```
Zhang Q., et al. (2025). GSIMEX: Generalized SIMEX with polynomial extrapolation. 
Code: https://github.com/QihuangZhang/GSIMEX
```

---

## Contact
Questions or issues? Please open a GitHub Issue or contact the maintainer.
