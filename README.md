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

# Pipeline Configuration

## Analysis Parameters

For specific parameters, file types, and other details needed for [regenie](https://rgcgithub.github.io/regenie/), [LDstore2](http://christianbenner.com/#), and [susieR](https://stephenslab.github.io/susieR/reference/susie_rss.html), please refer to their respective documentations. 

- `full_pipeline`: If set to `false`, only the `LDstore` and `SuSiE` components of the pipeline are executed.
- `regenie_outputs`: A text file (`.txt`) where each row corresponds to the `regenie` output file that should be used as input for `LDstore` and `SuSiE` when the previous parameter was set to `false`.

### Parameters for regenie step I

- `step1_input_format`: Select the input format for `regenie` step 1: `bed`, `pgen`, or `bgen`.
- `step1_bed`: When selecting `bed` as the input format, it assumes that the input folder contains the required files with the extensions `.bed`, `.bim`, and `.fam`.
- `step1_pgen`: When selecting `pgen` as the input format, it assumes that the input folder contains the required files with the extensions `.pgen`, `.pvar`, and `.psam`.
- `step1_bgen`: When selecting `bgen` as the input format, it assumes that the input folder contains the required files with the extensions `.bgen`, `.sample`, and `.bgi`.

### Parameters for regenie step I and II

- `phenotype_file`: Path to the phenotype file.
- `phenotype_list`: Path to a text file containing a comma-separated list of phenotypes.
- `covariate_file`: Path to the covariate file.
- `covariate_list`: A comma-separated list of covariates.
- `bsize`: Size of the genotype blocks. The default value is 4000.

## Parameters for regenie step II

- `samples_to_keep`: Path to a file specifying samples to keep.
- `minINFO`: Minimum imputation info score. The default value is 0.6.

## Parameters for regenie step II and LDStore2

All the following files for regenie step 2 must adhere to the format `chr*.bgen`/`chr*.sample`/`chr*.bgen.bgi`, where `*` represents the chromosome number.
- `step2_bgen`: Path to BGEN files for regenie step 2.
- `step2_sample`: Path to sample files for regenie step 2.
- `step2_bgi`: Path to BGI files for regenie step 2.

### Parameters for region filter

- `p_value`: p-value (default is 5e-08) for defining the genome-wide significant locus.
- `window_size`: Half of the size of the fine-mapping region. The default value is 1500000, which corresponds to a 3 Mb total window around a lead variant.
- `max_region_width`: Maximal fine-mapping region width (default is 6 Mb).
- `window_shrink_ratio`: Shrink ratio (default is 0.95) applied when regions exceed the maximal width (original regions are reduced recursively).
- `bgen_chr_has_zero`: Should be set to `T` (true) if chromosomes `1..9` are given as `01..09` in bgens.
- `remove_MHC`: If set to `T` (true), then the histocompatibility complex (MHC) region is excluded from the analysis.
- `GRCh`: Select the Right GRCh (choose either `37` or `38` for the GRCh).
- `MHC_start` and `MHC_end` : MHC Regions in GRCh37 and GRCh38 (specify the MHC region coordinates).

### Parameters for LDStore2

- `samples_ld_incl`: Path to the text file with samples that should be included.

### Parameters for susieR

- `n_covariates`: Number of covariates used in the analysis.
- `max_causal_SNPs`: Number of maximum causal variants (default is 10).

### Configuring output directory names

- `outdir`: Directory where results are stored (default is 'results').
- `prefix1`: Prefix for regenie results (default is 'regenie').
- `prefix2`: Prefix for fine-mapping results (default is 'finemapping').

## Execution configuration

- `excutor`: Process Executor
- `queue`: Process Queue

## Additional configuration

- Include the `base.config` file.

## Singularity Configuration

- Singularity is enabled.
- Auto-mounts are enabled.
- Cache directory: `"$baseDir/singularity_img/"`

## Helper Functions

- `check_max`: A function to check and validate max memory, time, and CPU parameters.
