# Vignette for using LOCATER for associaiton study

This repository contains the code used for the paper entitled Genealogy based trait association with LOCATER boosts power at loci with allelic heterogeneity. It also serves as a collection of vignettes for how to use LOCATER, a trait association procedure based on inferred local ancestries.

The steps from genotype data to association signal could be summarized in the following diagram:

![xinxin_Aug_23_2024 001](https://github.com/user-attachments/assets/13c23064-4cc6-4f6e-bc20-5f449f9de1a8)

For reference, the publication for LOCATER procedure is: [cite LOCATER]

The GitHub repository for the code of LOCATER's original paper is [refer to Ryan's LOCATEr paper github]

The original Github repository for LOCATER software is: [here](https://github.com/ryanchrist/locater)

## Dependencies

The procedure is based on `kalis` and `locater`, R packages developed by our group. For installation and vignettes of `kalis`, see [here](https://kalis.louisaslett.com/index.html).

Other dependencies include `dplyr`, `RSpectra`, `parallel`, `glmnet`, `tidyr`, `data.table`, and `Matrix`.

## Steps

In the following sections, we will describe the details of the LOCATER pipeline. The input files are as follows:
1. genotypes in the format of hap/sample/legend. The genotype files need to be coded so that the ancestral allele is REF allele, and phasing need to be performed.
2. (Optional) if phasing was done in segments, we need a table for the current segment boundaries of genotype files. The table need to have some way to query the genotype file (e.g. a "filename" column) and the core region start and end position (`core.start` and `core.end` column)
3. background covariate matrix, with columns as covariates
4. phenotype residuals, after accounting for phenotype-specific covariates
5. map files, the required columns are REF, ALT and map. We suggest using a table with columns chr, POS (genomic position), ID(any form of variant ID), REF, ALT, MAF, AC (allele count), AF, cM (genetic map in the unit of cM)

### Step 1: Tuning

Tuning steps correspond to scripts in `tuning` folder. 

#### Step 1.1: Run tuning

Run tuning corresponds to scripts `/tuning/submit.sh`, `/tuning/run_tuning`, and `/tuning/tuning_pars/METSIM1.R` See README inside  `tuning` folder for more information.

Input files of this step are: genotypes (the hap.gz file), segment boundary table with the ability to query the hap.gz file and the core region. 

Output files are `.rds` files containing results. The R object inside this file is a list containing the following elements:  "res", "smt.res" and diagnostic elements. 

`res` is a table containg tuning results from LOCATER, with `neglog10Ne` and `neglog10mu` as the two parameters to tune, and `rd`, `qform` and `signal` as the SD signal, QForm signal and combined signal, respectively. 

`smt.res` is a list, the length of this list is the number of parameter combinations tested.  We will query each of these elements during data collection to normalize LOCATER results. 

#### Step 1.2: Collect tuning data

Collect tuning data corresponds to scripts `/tuning/tuning_data_collect.R`

Input files of this step are: .rds` files containing results.

Output file is a  `.txt` file with normalized SD, QForm and combined signal. Other diagnostic columns contain derived allele count for the causal variant (DAC column), and the signal of SMT at target and causal variant (smt.target and smt.causal columns)

#### Step 1.3: Visualize tuning data

Visualize tuning data corresponds to scripts `/tuning/tuning_viz.R`. We suggest running this stem locally in RStudio.

Input files of this step are: `.txt` file from the last step.

Output files are HTML files containing interactive plots. By looking at the plots, we would also have a conclusion about the best parameters.

### Step 2: Preprocessing

Preprocessing steps correspond to scripts in the `preprocess` folder. 


#### Step 2.1: Define practical small segments

For efficiency, we suggest users of LOCATER separate the whole genome into smaller segments. In this specific study, during tuning we separated the whole genome into 155 segments, and I separated that into 4587 small segments. To ensure the running time for each segment is similar, we separate the genome based on the number of variants. 

Define practical small segments corresponds to scripts `/preprocess/1-define-segments.R`.

Input files of this step are: segment boundary table, target number of small segments

Output files is a table containing the boundaries of small segments. This table will contain columns with the start and end position of the core region for small segments, and the number of variant inside each small segment.

### Step 3: Rank matching and selection

Rank matching and selection step correspond to scripts in `rank_matching_selection` folder. 

#### Step 3.1: Generate rank matched phenotypes

Generate rank matched phenotypes corresponds to scripts `rank_matching_selection/1-generate-rank-matched-phenos.R`

Input files of this step are: phenotype residuals, after accounting for phenotype-specific covariates & number of rank matched versions for each phenotype. As mentioned in the script, we recommend setting a seed for reproducibility.

Output files are `.rds` files containing rank matched phenotypes. The R object inside this file would be a table with column number correspond to phentoype * version number, and row number correspond to sample size. 

#### Step 3.2: P-value calculation for rank matched phenotypes

Run tuning corresponds to scripts inside `rank_matching_selection/2-pval-rank-matched/` folder, See README inside  `rank_matching_selection` folder for more information.

Input files of this step are: .rds` files containing rank matched phenotypes, hap.gz files, table containing the boundaries of small segments, background covariate matrix, best HMM parameter setting.

Output files are `.rds` files containing p-value from LOCATER and its sub-tests, and the basic information about this segment. We will parse these files in the next step.

#### Step 3.3: rds to txt conversion

rds to txt conversion corresponds to scripts `/rank_matching_selection/3-rds-to-txt.R`. 

Input files of this step are: `.rds` files containing p-value from LOCATER and its sub-tests.

Output files are `.txt` files containing p-value from LOCATER and its sub-tests. The txt file looks like this:

```sh

locater.pos	locus.idx	thresh	max1var	old.sprigs	smt.noise	max.k	sw.thresh	eig.thresh	calc.obs.T	test.config	num.sprigs	k	exit.status	precise	obs.qform	obs.qform.T	smt	rd	qform	phenotype	num.layers	tot
100660841	4708	0.2	TRUE	FALSE	raw	512	0	0	TRUE	1	1715		0	TRUE	2384.13873736195	3.83509192902334	0.0142950705965422	0.291522825480577	0.147509569391811	eGFR_MDRD_combined1	1715	0.0148665949101459
100660841	4708	0.2	TRUE	FALSE	raw	512	0	0	TRUE	1	1715		0	TRUE	2402.27636009486	5.80228294940901	0.0296405375635567	0.174713259901607	0.2014313935507	eGFR_MDRD_combined2	1715	0.00637526569523721
100660841	4708	0.2	TRUE	FALSE	raw	512	0	0	TRUE	1	1715		0	TRUE	2399.36250921897	6.52074977312461	0.0202325851328748	0.216616778557003	0.191983104110816	eGFR_MDRD_combined3	1715	0.0068832040751616
100660841	4708	0.2	TRUE	FALSE	raw	512	0	0	TRUE	1	1715		0	TRUE	2429.82050989489	5.36201879833325	0.0108321781611223	0.523442750494169	0.306542692943568	eGFR_MDRD_combined4	1715	0.0668638754900115
100660841	4708	0.2	TRUE	FALSE	raw	512	0	0	TRUE	1	1715		0	TRUE	2440.32363063338	5.11069361905978	0.0239702916491983	0.314236858677153	0.354482835461624	eGFR_MDRD_combined5	1715	0.0272086773925919
100660841	4708	0.2	TRUE	FALSE	raw	512	0	0	TRUE	1	1715		0	TRUE	2441.21178034555	5.76999838439777	0.0199828714180033	0.141261994702723	0.358742522117343	eGFR_MDRD_combined6	1715	0.0272621550682296
100660841	4708	0.2	TRUE	FALSE	raw	512	0	0	TRUE	1	1715		0	TRUE	2357.12224216727	5.30993150169523	0.0175710795094786	0.138830141125491	0.0872553305419177	eGFR_MDRD_combined7	1715	0.00179456706086439
100660841	4708	0.2	TRUE	FALSE	raw	512	0	0	TRUE	1	1715		0	TRUE	2442.51941288333	6.04158521059452	0.0197504950384146	0.232864615132434	0.36507306022039	
...
```

The most important columns are `locater.pos` (genomic position),  `smt`	`rd` and `qform` (-log10(p-value) from SMT, SD and QForm, respectively) 

#### Step 3.4: Calculate generalized genomic control parameters & select the best rank matched version

Calculate generalized genomic control parameters & select the best rank matched version corresponds to scripts `/rank_matching_selection/4-calc-general-lambda-selection.R`. Note that this would only choose the best rank matched version if it exists. If there is no version that meet the requirement, the user will need to use the rank normalized version or choose other methods. 

Input files of this step are: `.txt` files from the last step.

Output files are `.rds` files containing all p-values from all versions of all phenotypes for LOCATER tests + slope and intercept table for all versions of all phenotypes + best rank matched version of each phenotypes + slope and intercept of SD and QForm for them.

Example output table:

```sh
qform_slope	qform_intercept	qform_rsquared	sd_slope	sd_intercept	sd_rsquared
eGFR_MDRD_combined1	0.508820670725126	-0.0326108773149501	0.997296962730409	0.795053687877007	0.476159092255134	0.996914918466741
eGFR_MDRD_combined2	0.686269541604983	-0.0338662266458304	0.994240749568966	0.832026767930068	0.169269291328179	0.986373188825943
eGFR_MDRD_combined3	0.611873055274072	-0.0660477857664275	0.996600438414666	0.86764950954921	0.00958478183543478	0.99681814098545
eGFR_MDRD_combined4	0.641369402472646	0.288364659612988	0.998235487645286	0.743312615397704	0.369609374215096	0.995957292121912
eGFR_MDRD_combined5	0.651910387257381	0.0841857270444568	0.99716062026905	0.970337470459158	0.102075448603857	0.985612695714009
eGFR_MDRD_combined6	0.862954269502971	0.131315337256212	0.997778329166161	0.974714797098509	0.0326949593364232	0.997007119309646
eGFR_MDRD_combined7	0.517148848880036	-0.137874043205288	0.996863896899126	0.991041598565586	0.0902682833851026	0.989864747359683
eGFR_MDRD_combined8	0.645229869883769	0.579357982064461	0.99838191959055	0.762283151538509	0.320563674658935	0.989685748175472
eGFR_MDRD_combined9	0.726160794594468	0.0228578576113133	0.996854451330703	0.880722242784847	0.0244290822625255	0.996713895091436
eGFR_MDRD_combined10	0.727687655043438	0.168011271638041	0.997854149015591	0.80042055492133	0.132836554256918	0.997251479450917
eGFR_MDRD_combined11	0.60969345979684	-0.0394305104755236	0.998527945145333	0.900108449797393	-0.0286895685995867	0.989203495690899
eGFR_MDRD_combined12	0.708456613533558	0.187668728978291	0.997465923347096	1.07757796231562	-0.122697584927048	0.993771522458657
eGFR_MDRD_combined13	0.672578165427807	-0.131471857279352	0.998272966537901	0.88338301073708	0.155363178448563	0.994284409436878
eGFR_MDRD_combined14	0.872166434118532	0.380512232570126	0.998228667327908	0.826385173476435	0.196864478306705	0.997307768853606
......
```

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
