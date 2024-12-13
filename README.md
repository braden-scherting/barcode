# barcode

Code to reproduce analyses within manuscript titled "Inferring latent structure in ecological communities via barcodes""

# Structure

- `SetupFunctions.R` supplies utility functions
- `ParallelFit.R` and `ParallelFitSpatial.R` 
- `progressFiles/` will be populated by files documenting sampling progress when the Gibbs sampler is run. The file `out0.txt` will be created and updated here when running a single chain, and files `out{X}.txt` will be created when running parallel chains $X\in1:n_{chains}$