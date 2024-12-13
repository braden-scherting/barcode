source("./SetupFunctions.R")

tmp <- loadData(noFactors=8)

runOneChainSpatial <- function(dims, dat, storeInit, workerSeed=sample(10:1e7,1), 
                        workerNo=0, niter=200, nthin=1, nburn=200, pred=0, blockSize=5){
  sourceCpp("./models/sampler_snapdragonSpatial.cpp")
  
  blockSize <- min(dims$L, blockSize)
  storeInit$Cstar <- matrix(0, ncol=blockSize, nrow=2^blockSize)
  for (i in 1:(2^blockSize - 1)){
    storeInit$Cstar[i+1,] <- binaryLogic::as.binary(i, n=blockSize)
  }
  
  return(sampler_snapdragonSpatial(niter, nthin, nburn, dims, dat, storeInit, 
                            workerSeed, x=workerNo, pred=pred))
}

runParallelSpatial <- function(inputs, nParallelChains=4, seedSeed=1292024, ...){
  
  registerDoFuture()
  plan(multisession, workers = nParallelChains)
  
  predTested <- ifelse(any(is.na(inputs$dat$Y)) && pred==TRUE, 1, 0)
  
  tmpFit <- foreach(i = 1:4, .options.future = list(seed = TRUE)) %dofuture% {
    runOneChainSpatial(dims=inputs$dims, dat=inputs$dat, storeInit=inputs$storeInit, 
                workerSeed=round(seedSeed * i / (i+1)), workerNo=i, pred=predTested, 
                ...)
  } 
  
  saveRDS(tmpFit, paste0("results/spatial", format(Sys.time(), "%H%M%a%d%b%Y"), ".rds"))
  return(NULL)
}

runParallelSpatial(tmp, niter=1500, nburn=70000, nthin=20, seedSeed = 8342024)