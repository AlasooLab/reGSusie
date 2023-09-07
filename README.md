# reGSusie Nextflow Pipeline
![GitHub](https://img.shields.io/badge/Version-1.0.0-blue.svg)
![GitHub](https://img.shields.io/badge/Nextflow-22.04.3-brightgreen.svg)

![pipeline_full](https://github.com/idarahu/reGSusie/assets/102286655/e61a86a3-773f-45d8-a89f-0b0748a45da4)


## Overview

The `reGSusie` Nextflow pipeline is designed for conducting Genome-Wide Association Studies (GWAS) using a combination of [regenie](https://rgcgithub.github.io/regenie/), [LDstore2](http://christianbenner.com/#), and [susieR](https://stephenslab.github.io/susieR/reference/susie_rss.html). This pipeline offers a comprehensive and customizable solution for processing genetic association data and fine-mapping causal variants.

### Key Features

- Conduct association testing with regenie.
- Compute in-sample LD matrix with LDstore2.
- Perform statistical fine-mapping with susieR.

## Prerequisites

Before running the pipeline, ensure you have the following dependencies installed:

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)
- [Singularity](https://sylabs.io/singularity/)

## Running the Pipeline

1. Clone the reGSusie GitHub repository and navigate to the directory:

   ```shell
   git clone https://github.com/idarahu/reGSusie.git
   cd reGSusie
2. Customize the nextflow.config file based on your analysis requirements (see Pipeline Configuration below).
3. To execute the pipeline, use the `submit_slurm.sh` script. This script is configured to submit a Slurm job for pipeline execution. Here's how to run the pipeline:
   ```bash
   sbatch submit_slurm.sh
