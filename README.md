# barcode

Code to reproduce analyses within manuscript titled *Joint species distribution modeling of abundance data through latent variable barcodes*.

# Structure


- `data/birds2.rds` Dataset; details below.
- `output/progressFiles/` contains blank text files to which MCMC progress updates are written
- `output/results`
  - This and that
- `R/` contains R scripts to reproduce analyses, simulations and figures. `barcode/R` is the expected working directory.
  - `R/setup_functions.R` supplies utility functions and is called from other scripts. 
  - Fit models using `R/analysis.R` in an R interactive console. 
  - Generate figures from Section 4 using `R/plotting_functions.R`.
  - Compute statistical summaries referenced in Section 4 using `R/summaries.R`.
  - `R/simulations/` contains scripts for performing simulations and plotting results described in the Supplementary Materials. 
- `src/` contains `.cpp` implementations of the model. 
- `data/` contains the data.
- `output/` contains two folders; `output/progressFiles` holds text files that are updated during sampling for monitoring progress, and `output/results` hold fitted models and simulation results. 

# Data Dictionary 

The `birds2.rds` dataset is a list with the following elements:

- `Y` is a 2826-by-132 matrix of bird abundances. Row and column names identify the sample ID and species name, respectively.
- `XData` is a 2826-by-21 matrix of covariates
  - `Mixed_Forest`
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