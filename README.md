# barcode

Code to reproduce analyses within manuscript titled ``Inferring latent structure in ecological communities via barcodes''

# Structure

\begin{itemize}
  \itm \verb|SetupFunctions.R| supplies utility functions
  \item \verb|ParallelFit.R| and \verb|ParallelFitSpatial.R| 
  \item \verb|progressFiles/| document sampling progress when the Gibbs sampler is run. The file \verb|out0.txt| will be created and updated here when running a single chain, and files \verb|out{X}.txt| will be created when running parallel chains $X\in1:n_{chains}$
\end{itemize}