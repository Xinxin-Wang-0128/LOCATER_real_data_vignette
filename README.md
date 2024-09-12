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
4. phenotype residuals, after accounting for phenotype-specific covariates
5. map files, the required columns are REF, ALT and map

### Step 1: Tuning

Tuning steps correspond to scripts in `tuning` folder. 

#### Step 1.1: Run tuning

Run tuning corresponds to scripts `/tuning/submit.sh`, `/tuning/run_tuning`, and `/tuning/tuning_pars/METSIM1.R` See README inside  `tuning` folder for more information.

Input files of this step are: genotypes, segment boundary table. 

Output files are `.rds` files containing results. [add detail of this file format]

#### Step 1.2: Collect tuning data

Collect tuning data corresponds to scripts `/tuning/tuning_data_collect.R`

Input files of this step are: .rds` files containing results.

Output files are `.txt` files [add detail of this file format]

#### Step 1.3: Visualize tuning data

Visualize tuning data corresponds to scripts `/tuning/tuning_viz.R`. We suggest running this stem locally in RStudio.

Input files of this step are: `.txt` files from the last step

Output files are HTML files containing interactive plots. [add detail of this file format]. We would also have a conclusion about the best parameters.

### Step 2: Preprocessing

Preprocessing steps correspond to scripts in the `preprocess` folder. 


#### Step 2.1: Define practical small segments

For efficiency, we suggest users of LOCATER separate the whole genome into smaller segments. In this case, after tuning we have 155 segments, and I separated that into 4587 segments. To ensure the running time for each segment is similar, we separate the genome based on the number of variants. 

Define practical small segments corresponds to scripts `/preprocess/1-define-segments.R`.

Input files of this step are: segment boundary table, target number of small segments

Output files is a table containing the boundaries of small segments. [add detail of this file format]

### Step 3: Rank matching and selection

Rank matching and selection step correspond to scripts in `rank_matching_selection` folder. 

#### Step 3.1: Generate rank matched phenotypes

Generate rank matched phenotypes corresponds to scripts `rank_matching_selection/1-generate-rank-matched-phenos.R`

Input files of this step are: phenotype residuals, after accounting for phenotype-specific covariates & number of rank matched versions for each phenotype 

Output files are `.rds` files containing rank matched phenotypes. The R object inside this file would be a table with column number correspond to phentoype * version number, and row number correspond to sample size. [add detail of this file format]

#### Step 3.2: P-value calculation for rank matched phenotypes

Run tuning corresponds to scripts inside `rank_matching_selection/2-pval-rank-matched/` folder, See README inside  `rank_matching_selection` folder for more information.

Input files of this step are: .rds` files containing rank matched phenotypes, hap files, table containing the boundaries of small segments, background covariate matrix, best HMM parameter setting.

Output files are `.rds` files containing p-value from LOCATER and its sub-tests. [add detail of this file format]

#### Step 3.3: rds to txt conversion

rds to txt conversion corresponds to scripts `/rank_matching_selection/3-rds-to-txt.R`. 

Input files of this step are: `.rds` files containing p-value from LOCATER and its sub-tests.

Output files are `.txt` files containing p-value from LOCATER and its sub-tests. [add detail of this file format].

#### Step 3.4: Calculate generalized genomic control parameters & select the best rank matched version

Calculate generalized genomic control parameters & select the best rank matched version corresponds to scripts `/rank_matching_selection/4-calc-general-lambda-selection.R`. Note that this would only choose the best rank matched version if it exists. If there is no version that meet the requirement, the user will need to use the rank normalized version or choose other methods. 

Input files of this step are: `.txt` files from the last step

Output files are `.rds` files containing all p-values from all versions of all phenotypes for LOCATER tests + slope and intercept table for all versions of all phenotypes + best rank matched version of each phenotypes + slope and intercept of SD and QForm for them [add detail of this file format]. 

### Step 4: Whole genome screening

#### Step 4.1 Screening

Screening correspond to scripts in `whole_genome_screening/2-screening` folder. 

Input files of this step are: robust phenotypes (best rank matched version + rank normalized for phenotypes that cannot pick the best rank matched version); hap files, table containing the boundaries of small segments, background covariate matrix, best HMM parameter setting.

Output files are `.rds` files containing p-value from LOCATER and SMT. [add detail of this file format]

#### Step 4.2: rds to txt conversion

rds to txt conversion corresponds to scripts `/whole_genome_screening/3-rds-to-txt.R`. 

Input files of this step are: `.rds` files containing p-value from LOCATER and SMT.

Output files are `.txt` files containing p-value from LOCATER and SMT. [add detail of this file format].

#### Step 4.3 Collect screening data and visualize interesting associations 

Collect screening data and find interesting associations correspond to script `whole_genome_screening/4-locater-wg-viz.R`. 

Input files of this step are: `.txt` files containing p-value from LOCATER and SMT and slope and intercept table of SD and QForm for robust phenotypes.

Output files are `.rds` files containing whole chromosome p-value from LOCATER and SMT and `.rds` files containing local p-values for interesting associaitons. [add detail of this file format]

We would also create PDF files containing interesting association.

#### Step 4.4 Collect information for interesting associations 

Collect information for interesting associations correspond to script `whole_genome_screening/5-find-interesting-assoc.R`. 

Input files are `.rds` files containing local p-values for interesting associations.

Output file is a table that contains basic information for all interesting associations

#### Step 4.5 Generate candidate loci 

Generate candidate loci correspond to script `whole_genome_screening/6-generate-candidate-loci.R`. 

Input file is the table that contains basic information for all interesting associations. [add detail of this file format]

Output file is a table for basic information for candidate loci.

### Step 5: Candidate loci investigation

#### Step 5.1 Exact p-value calculation

Screening correspond to scripts in `followup/3-putative-investigation` folder. Note that for all candidate loci we would want to set the same seed in R, and use this seed for all later investigations. 

Input files of this step are: robust phenotypes (best rank matched version + rank normalized for phenotypes that cannot pick the best rank matched version); hap files, background covariate matrix, best HMM parameter setting and basic information for candidate loci (chr, start and end)

Output files are `.rds` files containing p-value from LOCATER and SMT. [add detail of this file format] Users could use the same rds to txt conversion and data collection method described in step 4.2 and 4.3.

#### Step 5.2 Sub-test evaluation

This step is to visualize which sub-test contributed to LOCATER signal locally for the loci or specifically for the lead marker, and to what extend the adjustments and standardization changed the p-value. This would base on a `.rds` file that contain the p-value table for SMT and LOCATER for interested loci and a interested phenotype. 

The script is `followup/sub-test-evaluation.R` 
Output files are PDF files with Manhattan plots and bar charts.

#### Step 5.3: Sprig object generation

This is only useful when SD is the sub-test that boosted the LOCATER signal.

Sprig object generation corresponds to scripts `/followup/sprig_object_generation`. 

Input files of this step are: apart from the input files in step 5.1, we would also need the position of interested variant.

Output files are various `.rds` files containing .... [add detail of this file format] 

#### Step 5.4: Sprig object visualization

This is only useful when SD is the sub-test that boosted the LOCATER signal.

Sprig object visualization corresponds to scripts `/followup/sprig_object_viz.R` and `/followup/hap_matrix_viz.R`

Input files of this step are: output files from step 5.4

Output files are PDF files containing plots. 

#### Step 5.5: Residual analysis [may not want to add this]

This is only useful when SD is the sub-test that boosted the LOCATER signal.

Residual analysis corresponds to scripts `/followup/???`. [add this]

Input files of this step are: apart from the input files in step 5.1, we would also need the position of interested variant.[really??]

Output files are various `.rds` files containing .... [add detail of this file format] 
