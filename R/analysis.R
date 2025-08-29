source("setup_functions.R")

dat5 <- loadData(noFactors=5, initialization = T)
runParallel(dat5, niter=2500, nburn=25000, nthin=10, seedSeed = 320104,
            nParallelChains=4, prefix="rank5_")

dat6 <- loadData(noFactors=6, initialization = T)
runParallel(dat6, niter=2500, nburn=25000, nthin=10, seedSeed = 424257,
            nParallelChains=4, prefix="rank6_")

dat7 <- loadData(noFactors=7, initialization = T)
runParallel(dat7, niter=2500, nburn=25000, nthin=10, seedSeed = 947390,
            nParallelChains=4, prefix="rank7_")

dat8 <- loadData(noFactors=8, initialization = T)
runParallel(dat8, niter=2500, nburn=25000, nthin=10, seedSeed = 509262,
            nParallelChains=4, prefix="rank8_")

dat10 <- loadData(noFactors=10, initialization = T)
runParallel(dat10, niter=2500, nburn=25000, nthin=10, seedSeed = 625114,
            nParallelChains=4, prefix="rank10_")

dat20 <- loadData(noFactors=20, initialization = T)
runParallel(dat20, niter=2500, nburn=25000, nthin=10, seedSeed = 161152,
            nParallelChains=4, prefix="rank20_")
