# GWAS_Nextflow_pipeline

![reGSusie](https://github.com/idarahu/reGSusie/assets/102286655/b838d762-29ba-473d-9063-a52a5f4d82ae)

![reGSusie](https://github.com/idarahu/reGSusie/assets/102286655/148b0606-f24c-4c15-871f-8cd7ee9716e3)# reGSusie Nextflow Pipeline

![GitHub](https://img.shields.io/badge/Version-1.0.0-blue.svg)
![GitHub](https://img.shields.io/badge/Nextflow-22.04.3-brightgreen.svg)

## Overview

The `reGSusie` Nextflow pipeline is designed for conducting Genome-Wide Association Studies (GWAS) using a combination of REGENIE, LDstore2, and SuSiE. This pipeline offers a comprehensive and customizable solution for processing genetic association data and fine-mapping causal variants.

### Key Features

- Conduct GWAS analysis using REGENIE.
- Perform fine-mapping with LDstore2.
- Apply SuSiE for Bayesian sparse regression.
- Flexible configuration options for different input data formats and analysis scenarios.

## Prerequisites

Before running the pipeline, ensure you have the following dependencies installed:

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)
- [Singularity](https://sylabs.io/singularity/)

## How to Run

To execute the `reGSusie` pipeline, follow these steps:

1. Clone the repository:

   ```shell
   git clone https://github.com/idarahu/reGSusie.git
   cd reGSusie
