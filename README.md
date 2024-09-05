# Vignette for using LOCATER for associaiton study

This repository contains the code used for the paper entitled Genealogy based trait association with LOCATER boosts power at loci with allelic heterogeneity. It also serves as a collection of vignettes for how to use LOCATER, a trait association procedure based on inferred local ancestries.

The steps from genotype data to association signal could be summarized in the following diagram:

![xinxin_Aug_23_2024 001](https://github.com/user-attachments/assets/13c23064-4cc6-4f6e-bc20-5f449f9de1a8)

For reference, the publication for LOCATER procedure is: [cite LOCATER]

The GitHub repository for the code of LOCATER's original paper is: [refer to Ryan's GitHub repo]

The original Github repository for LOCATER software is: [refer to LOCATER github]

## Dependencies

The procedure is based on `kalis` and `locater`, R packages developed by our group. For installation and vignettes of `kalis`, see [here](https://kalis.louisaslett.com/index.html).

Other dependencies include `dplyr`, `RSpectra`, `parallel`, `glmnet`, `tidyr`, `data.table`, and `Matrix`.

## Steps

In the following sections, we will describe the details of the LOCATER pipeline. The input files are as follows:
1. genotypes in the format of hap/sample/legend. The genotype files need to be coded so that the ancestral allele is REF allele, and phasing need to be performed.
2. (Optional) if phasing was done in segments, we need a table for the current segment boundaries of genotype files.
3. background covariate matrix, with columns as covariates
4. phenotype residuals
5. map files, the required columns are REF, ALT and map

### Step 1: Tuning

Tuning steps correspond to scripts in `tuning` folder. 

#### Step 1.1: Run tuning

Run tuning corresponds to scripts `/tuning/submit.sh`, `/tuning/run_tuning`, and `/tuning/tuning_pars/METSIM1.R` See README inside  `tuning` folder for more information.

Input files of this step are: genotypes, segment boundary table. 

Output files are `.rds` files containing results. [add detail of this file format]

#### Step 1.2: Collect tuning data

Run tuning corresponds to scripts `/tuning/tuning_data_collect.R`

Input files of this step are: .rds` files containing results.

Output files are `.txt` files [add detail of this file format]

#### Step 1.3: Visualize tuning data

Run tuning corresponds to scripts `/tuning/tuning_viz.R`. We suggest running this stem locally in RStudio.

Input files of this step are: `.txt` files from the last step

Output files are HTML files containing interactive plots. [add detail of this file format]. We would also have a conclusion about the best parameters.

### Step 2: Preprocessing

Preprocessing steps correspond to scripts in the `preprocess` folder. 


#### Step 2.1: Define practical small segments

For efficiency, we suggest users of LOCATER separate the whole genome into smaller segments. In this case, after tuning we have 155 segments, and I separated that into 4587 segments. To ensure the running time for each segment is similar, we separate the genome based on the number of variants. 

Define practical small segments corresponds to scripts `/preprocess/1-define-segments.R`.

Input files of this step are: segment boundary table, target number of small segments

Output files is a table containing the boundaries of small segments. [add detail of this file format]





