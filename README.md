# barcode

Code to reproduce analyses within manuscript titled *Joint species distribution modeling of abundance data through latent variable barcodes*.

# Structure


- `data/birds2.rds` Dataset; details below.
- `output/progressFiles/` contains blank text files to which MCMC progress updates are written
- `output/results/`
  - `recoverC.rds` contains fitted model object for reproducing Figure A.4.
  - `recoverS.rds` contains fitted model object for reproducing Figure A.4.
  - `BigBeta2.rds` contains parameter estimates for reproducing Figure A.5.
  - `OOPSEpreds.rds`contains cross validation predictions to reproduce results presented under the heading ``Explaining and Predicting Observed Community Distribution''
  - `foldRMSE.rds` contains summaries based on `OOPSEpreds.rds`
- `R/` contains R scripts to reproduce analyses, simulations and figures. `barcode/R` is the expected working directory.
  - `R/setup_functions.R` is called from other scripts and supplies utility functions for loading data in the expected format, model fitting, and compiling results. 
  - Fit models using `R/analysis.R` in an R interactive console. The provided model fit calls produce results for $L\in\{5,6,7,8,10,20\}$.
  
  - Generate figures from Section 4 using `R/plotting_functions.R`. Manuscript figures can be replicated using this script based on the provided contents of `output/results`.
  - `R/summaries.R` included code to replicate numerical summaries based on contents of `output/results`.
  - `R/simulations/` contains scripts for performing simulations and plotting results described in the Supplementary Materials. 
    - `CSSparsityRecovery2.R` reproduces Figure A.4.
    - `BetaRecovery2.R` reproduces Figure A.5
    - `OutOfSamplePrediction2.R` reproduces results presented under the heading ``Explaining and Predicting Observed Community Distribution''
- `src/barcodeSampler.cpp` contains model fitting code and is called by other functions in `setup_functions.R`

# Data Dictionary 

The `birds2.rds` data set is a list with the following elements:

- `Y` is a 2826-by-132 matrix of bird abundances. Row and column names identify the sample ID and species name, respectively.
- `XData` is a 2826-by-21 matrix of covariates
  - `Mixed_Forest`: Fraction of area surrounding the site classified as mixed forest
  - `Deciduous_Forest`: Fraction of area surrounding the site classified as deciduous forest
  - `Shrubs`: Fraction of area surrounding the site classified as shrubs
  - `Grasslands_Wetlands`: Fraction of area surrounding the site classified as grasslands OR wetlands
  - `Acricultural_Land`: Fraction of area surrounding the site classified as agricultural land
  - `Barren`: Fraction of area surrounding the site classified as barren
  - `Urban`: Fraction of area surrounding the site classified as urban
  - `Water_Bodies`: Fraction of area surrounding the site classified as water
  - `Coastal`: Fraction of area surrounding the site classified as coastal
  - `Stand_Age`: Age of the surrounding forest in years
  - `Pine_Volume`: Forest inventory metric: pine volume
  - `Spruce_Volume`: Forest inventory metric: spruce volume
  - `Birch_Volume`: Forest inventory metric: birch volume
  - `Other_Deciduous_Volume`: Forest inventory metric: other deciduous volume
  - `VMI_PC1`: First principal component of vegetation indices
  - `duration`: Survey duration (log)
  - `linelength`: Line transect length (log)
  - `polyAprMay1`: Annual periodic function peaking between April and May (sampling times to not cover the whole year)
  - `polyDecFeb1`: Annual periodic function peaking in January
  - `polyJunJul1`: Annual periodic function peaking between June and July
- `des01` is a 2826-by-555 binary matrix with one-hot rows indicating which site each sample comes from
- `D` is a 555-by-555 pairwise site-to-site distance matrix 
- `XFormula` is an `R` formula object describing the latent regression specification
- `long.lat` is a 555-by-2 matrix of site coordinates
- `year` is a 2826-vector stating the number of years prior to 2016 in which each sample was collected

# Scope

The provided workflow reproduces:

- Any numbers provided in the text of the paper
- The computational methods presented in the paper 
- All tables and figures in the paper

# Workflow

- The workflow is available at https://github.com/braden-scherting/barcode as wrapper scripts with further instructions provided below.

# Instructions

- To replicate main analysis in manuscript Section 4, execute `analysis.R` from the `barcode/R` working directory in an `R` console. MCMC samples will be saved at: `ouput/results/`.
- Figures 1–5 and A.1–A.3 use the saved objects and can be recreated using `R/plotting_functions.R`; update the name of the saved samples in each script as needed if model is re-fit.
- To replicate results from each simulation setting, execute the corresponding R script from either console or IDE:
  - Execute `CSSparsityRecovery2.R` to replicate sparse latent variable recovery; Figure A.4. Figure creation within.
  - Execute `OutOfSamplePrediction.R` to replicate out-of-sample predictions on case study data (bird abundances). Figure creation within.
- Running `R/summaries.R` gives the various numerical values referenced throughout the results section.
- The expected run-time is $>8$ hours