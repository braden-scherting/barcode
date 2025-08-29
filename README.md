# barcode

Code to reproduce analyses within manuscript titled *Joint species distribution modeling of abundance data through latent variable barcodes*.

# Structure

- `R/` contains R scripts to reproduce analyses, simulations and figures. `barcode/R` is the expected working directory.
  - `R/setup_functions.R` supplies utility functions and is called from other scripts. 
  - Fit models using `R/analysis.R` in an R interactive console. 
  - Generate figures from Section 4 using `R/plotting_functions.R`.
  - Compute statistical summaries referenced in Section 4 using `R/summaries.R`.
  - `R/simulations/` contains scripts for performing simulations and plotting results described in the Supplementary Materials. 
- `src/` contains `.cpp` implementations of the model. 
- `data/` contains the data.
- `output/` contains two folders; `output/progressFiles` holds text files that are updated during sampling for monitoring progress, and `output/results` hold fitted models and simulation results. 