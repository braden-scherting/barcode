source("SetupFunctions.R")

tmp <- loadData(noFactors=6)

runParallel(tmp, niter=1500, nburn=25000, nthin=50, seedSeed = 1092024)