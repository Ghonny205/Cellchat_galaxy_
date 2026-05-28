# CellChat Galaxy Wrapper

[![Galaxy](https://img.shields.io/badge/Galaxy-Tool-blue.svg)](https://galaxyproject.org/)
[![Docker](https://img.shields.io/badge/Docker-Containerized-blueviolet.svg)](https://www.docker.com/)

This repository contains a Galaxy tool wrapper for the [CellChat](https://github.com/jinworks/CellChat) R package. It allows users to run cell-cell communication (CCC) analysis on scRNA-seq data directly through the Galaxy web interface, without needing to write any R code.

## Main features & fixes

If you've tried installing CellChat recently, you probably ran into some dependency issues. This wrapper handles a few of those headaches automatically:

* **Seurat v4/v5 compatibility:** The R script uses a `tryCatch` block to dynamically extract matrices whether your `.rds` uses Seurat v4 (`slots`) or the new Seurat v5 (`layers`). You don't need to reformat your data beforehand.
* **NMF package compilation:** The provided Dockerfile uses Mamba but compiles the `NMF` dependency directly from CRAN using base system compilers. This bypasses the common Fortran/C library conflicts on newer architectures.
* **Simple GUI:** Exposes the basic CellChat pipeline (database selection, metadata grouping, and network inference) through a standard Galaxy XML interface.

## Repository structure

* `docker/Dockerfile`: Base image setup (Miniforge3, R, Mamba, and CellChat dependencies).
* `scripts/cellchat_wrapper.R`: The R script that parses Galaxy inputs, runs the CellChat pipeline, and exports the results.
* `tool/cellchat.xml`: The Galaxy tool definition file.

## Testing locally with Planemo

If you want to test or modify the tool locally, you can use [Planemo](https://planemo.readthedocs.io/en/latest/).

**Prerequisites:**
* Docker installed and running.
* Planemo (`pip install planemo`).

**Steps:**
1. Clone this repo:
```bash
   git clone [https://github.com/tu_usuario/cellchat_galaxy.git](https://github.com/tu_usuario/cellchat_galaxy.git)
   cd cellchat_galaxy/tool
