# Vignette for using LOCATER for associaiton study

This repository contains the code used for the paper entitled Genealogy based trait association with LOCATER boosts power at loci with allelic heterogeneity. It also serves as a collection of vignettes for how to use LOCATER, a trait association procedure based on inferred local ancestries.

The steps from genotype data to association signal could be summarized in the following diagram:

![xinxin_Aug_23_2024 001](https://github.com/user-attachments/assets/13c23064-4cc6-4f6e-bc20-5f449f9de1a8)

For reference, the publication for LOCATER procedure is: [cite LOCATER]. Please cite this paper if you use LOCATER in your study.

The GitHub repository for the code of LOCATER's original paper is [refer to Ryan's LOCATEr paper github]

The original Github repository for LOCATER software is: [here](https://github.com/ryanchrist/locater)

## Dependencies

The procedure is based on `kalis` and `locater`, R packages developed by our group. 

For installation and vignettes of `kalis`, see [here](https://kalis.louisaslett.com/index.html). 

For installation of LOCATER, see [here](https://github.com/ryanchrist/locater) [make sure this has installation instructions]

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

Input files of this step are: `.rds` files containing results.

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
Based on the best rank matched version of available phenotypes and the rank normalized phenotypes for ones that don't have the best rank matched version, users could construct a robust phenotype table, with columns being all phenotypes interested and the rows being individuals. 

### Step 4: Whole genome screening

#### Step 4.1 Screening

Screening correspond to scripts in `whole_genome_screening/2-screening` folder. 

Input files of this step are: robust phenotypes (best rank matched version + rank normalized for phenotypes that cannot pick the best rank matched version); hap files, table containing the boundaries of small segments, background covariate matrix, best HMM parameter setting.

Output files are `.rds` files containing p-value from LOCATER and SMT, and the basic information about this segment. We will parse these files in the next step.

#### Step 4.2: rds to txt conversion

rds to txt conversion corresponds to scripts `/whole_genome_screening/3-rds-to-txt.R`. 

Input files of this step are: `.rds` files containing p-value from LOCATER and SMT.

Output files are `.txt` files containing p-value from LOCATER and SMT.

For the `.txt` file of LOCATER, it looks like the following table:

```sh

locater.pos  locus.idx  thresh  max1var  old.sprigs  smt.noise  max.k  sw.thresh  eig.thresh  cs.approx  test.config  num.sprigs  k  exit.status  precise  prop.var  var.ratio          smt                rd                  qform               phenotype              num.layers tot
119387217    4928       0.2     TRUE     FALSE       raw        0      6          6.5         FALSE      1            1884        0  0            FALSE    0         NA                 0.490663367251888  0.493882584265179   0.0347904622373371  eGFR_MDRD_combined92   1884       0.276766447542598   
119387217    4928       0.2     TRUE     FALSE       raw        0      6          6.5         FALSE      1            1884        0  0            FALSE    0         NA                 0.549478227289189  0.612385643245848   0.0257749225565219  S_krea_combined94      1884       0.360707138613825   
119387217    4928       0.2     TRUE     FALSE       raw        0      6          6.5         FALSE      1            1884        0  0            FALSE    0         NA                 0.299944411579717  0.0431535961619506  0.109091640585671   S_totalc_combined52    1884       0.137696420561351   
119387217    4928       0.2     TRUE     FALSE       raw        0      6          6.5         FALSE      1            1884        0  0            FALSE    0         NA                 0.208427898040653  0.190727845883129   0.122268869273039   ln_S_tottg_combined35  1884       0.0797823593041662  
119387217    4928       0.2     TRUE     FALSE       raw        0      6          6.5         FALSE      1            1884        0  0            FALSE    0         NA                 0.569059802960148  0.593313671106326   0.348724894638233   ln_P_CRP_combined      1884       0.361006369922754   
119387217    4928       0.2     TRUE     FALSE       raw        0      6          6.5         FALSE      1            1884        0  0            FALSE    0         NA                 0.334549337965149  0.83504448802095    0.537047076462959   height_combined34      1884       0.367147946025518   
119387217    4928       0.2     TRUE     FALSE       raw        0      6          6.5         FALSE      1            1884        0  0            FALSE    0         NA                 0.386968486772063  0.142858306516035   0.0937612683956253  bmi_combined87         1884       0.198765880725852   
119387217    4928       0.2     TRUE     FALSE       raw        0      6          6.5         FALSE      1            1884        0  0            FALSE    0         NA                 0.461611094048419  0.235232368396069   0.41876800544341    weight_combined58      1884       0.254445621096923   
119387217    4928       0.2     TRUE     FALSE       raw        0      6          6.5         FALSE      1            1884        0  0            FALSE    0         NA                 0.150978716652847  0.924513318776337   0.413816623200519   systbp_combined3       1884       0.310958241456534   

...
```
and for the ones of SMT, it looks like:

```sh

smt.pos    map               MAC   smt.p               phenotype
119386926  121.050295674509  419   0.190122908513672   eGFR_MDRD_combined92
119386995  121.050357742536  2510  0.0968067455618945  eGFR_MDRD_combined92
119387217  121.050557130045  12    0.490663367251888   eGFR_MDRD_combined92
119387284  121.050596553     2507  0.0802841672351469  eGFR_MDRD_combined92
119387314  121.050610952984  48    0.636572077027618   eGFR_MDRD_combined92
119387440  121.050671401     6524  0.0281836966411413  eGFR_MDRD_combined92
119387461  121.050680836489  6526  0.0138721922906406  eGFR_MDRD_combined92
119387575  121.050730482235  6     0.31888977384413    eGFR_MDRD_combined92
119387870  121.050858938649  205   0.304594375162493   eGFR_MDRD_combined92
...

```

#### Step 4.3 Collect screening data and visualize interesting associations 

Collect screening data and find interesting associations correspond to script `whole_genome_screening/4-locater-wg-viz.R`. 

Input files of this step are: `.txt` files containing p-value from LOCATER and SMT and slope and intercept table of SD and QForm for robust phenotypes.

Output files are `.rds` files containing whole chromosome p-value from LOCATER and SMT and `.rds` files containing local p-values for interesting associations. 

The whole chromosome file contains a list where all results for each phenotype is an element of this list. The length of this list is the number of phenotypes, and each element contains results from all variants in this chromosome for this specific phenotype. 

The local file contains a list where the element `smt.in.plot` and `locater.in.plot` each contain a table that has the same format as the `.txt` file mentioned in Step 4.2, but only for the interested phenotype and region.

We would also create PDF files containing Manhattan plots of interesting association.

#### Step 4.4 Collect information for interesting associations 

Collect information for interesting associations correspond to script `whole_genome_screening/5-find-interesting-assoc.R`. 

Input files are `.rds` files containing local p-values for interesting associations from the last step.

Output file is a table that contains basic information for all interesting associations.

The output table should look like this:
```sh
chr  phenotype        start      end        locater_tot       smt               max.locater.pos  max.smt.pos  sig.locater
1    Alb_combined2    107211635  108411635  7.10109803447498  7.58941467622041  107811635        107811635    FALSE
1    Alb_combined2    1767871    2967871    6.34402427105624  4.49488204299527  2367871          2367871      TRUE
1    Alb_combined2    239448945  240648945  6.57465405866869  7.05650148421113  240048945        240048945    FALSE
1    bmi_combined87   155303523  156538032  6.65523524921212  7.13807289815438  155911161        155911161    FALSE
1    bmi_combined87   237356129  238562928  6.6284636818973   7.11097234795622  237956319        237956319    FALSE
1    bmi_combined87   246943495  248362310  6.96314299115726  6.21153892147095  247762310        247543495    TRUE
1    FAw3_combined51  48622192   49822192   6.30437131487927  6.78289736609534  49222192         49222192     FALSE
1    FAw3_combined51  56352598   57552598   6.4397314144312   6.91992084064918  56952598         56952598     FALSE
1    FAw6_combined    61835572   63312025   8.66664420286488  9.17419908667314  62499950         62499950     FALSE
...
```

#### Step 4.5 Generate candidate loci 

Generate candidate loci correspond to script `whole_genome_screening/6-generate-candidate-loci.R`. 

Input file is the table that contains basic information for all interesting associations from the last step.

Output file is a table for basic information for candidate loci.

The output table will look like this:

```sh
chr  start     end
1    1767871   2967871
1    3828771   5028771
1    16265079  17465079
1    28663692  29863692
1    41125005  42325005
1    46503256  47703256
1    48622192  49822192
1    53581044  57130021
1    56352598  57552598
...

```

### Step 5: Candidate loci investigation

#### Step 5.1 Exact p-value calculation

Screening correspond to scripts in `followup/3-putative-investigation` folder. Note that for all candidate loci we would want to set the same seed in R, and use this seed for all later investigations. 

Input files of this step are: robust phenotypes (best rank matched version + rank normalized for phenotypes that cannot pick the best rank matched version); hap files, background covariate matrix, best HMM parameter setting and basic information for candidate loci (chr, start and end) from step 4.5.

Output files are `.rds` files containing p-value from LOCATER and SMT for this interested region. Users could use the same rds to txt conversion and data collection method described in step 4.2 and 4.3.

#### Step 5.2 Sub-test evaluation

This step is to visualize which sub-test contributed to LOCATER signal locally for the loci or specifically for the lead marker, and to what extend the adjustments and standardization changed the p-value. This would base on a `.rds` file that contain the p-value table for SMT and LOCATER for interested loci and an interested phenotype. 

The script is `followup/sub-test-evaluation.R`.

Output files are PDF files with Manhattan plots and bar charts.

Example plots:

![bar_peak_ln_M_HDL_TG_combined83_49038347-50253146](https://github.com/user-attachments/assets/45b076bb-6f30-48a0-bc8e-255739f031c7)

![split_Manhattan_ln_M_HDL_TG_combined83_49038347-50253146_locaters](https://github.com/user-attachments/assets/d04554a8-8179-46b7-9901-26ce5f959db6)

#### Step 5.3: Sprig object generation

This is only useful when SD is the sub-test that boosted the LOCATER signal.

Sprig object generation corresponds to scripts `/followup/sprig_object_generation`. 

Input files of this step are: addition to the input files in step 5.1, we would also need the position of interested variant.

Output files are various `.rds` files containing information for many aspects of sprig object.

`spigs` is a object that include sprig assignment information and number of sprigs.

`res.inside` is a list that include detailed SD results.  `res.inside$u` is a matrix containing p-values from all sprigs called at this variants with all phenotypes; `res.inside$p_value` is an array with p-values from SD test for all phenotypes at this variant.  

`g` is the genotype vector at this interested variant. `A`: background covariate matrix. `y`: phenotypes.


#### Step 5.4: Sprig object visualization

This is only useful when SD is the sub-test that boosted the LOCATER signal.

Sprig object visualization corresponds to scripts `/followup/sprig_object_viz.R` and `/followup/hap_matrix_viz.R`

Input files of this step are: output files from step 5.4

Output files are PDF files containing plots: 

<img width="678" alt="Screenshot 2024-09-16 at 3 41 59 PM" src="https://github.com/user-attachments/assets/3ea0f23c-1116-42e8-b02b-224483bd9904">
<img width="746" alt="Screenshot 2024-09-16 at 3 42 39 PM" src="https://github.com/user-attachments/assets/00ff9b2a-d500-4227-9bd3-c4b59b008785">

<img width="1329" alt="Screenshot 2024-09-16 at 3 42 58 PM" src="https://github.com/user-attachments/assets/9f138b6e-851b-4155-977d-30f3137ac959">

#### Comments for associations where QForm boosted LOCATER signal

Users could choose their preferred way to do conditional analysis (condition on the lead marker of SMT iteratively) and use preferred way to visualize. In this publication, we used in-house scripts to plot LocusZoom plots.

#### Step 5.5: Residual analysis [may not want to add this]

This is only useful when SD is the sub-test that boosted the LOCATER signal.

Residual analysis corresponds to scripts `/followup/???`. [add this]

Input files of this step are: apart from the input files in step 5.1, we would also need the position of interested variant.[really??]

Output files are various `.rds` files containing .... [add detail of this file format] 
