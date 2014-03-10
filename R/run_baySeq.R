run_baySeq <- function(counts, conds, cutoff, n, runID) {
    cl <- makeCluster(2, "SOCK")
    # Preparing variables for baySeq run
    NR1 = length(conds[conds == "N"])
    NR2 = length(conds[conds == "T"])
    groups <- list(NDE = c(rep(1, NR1 + NR2)), DE = c(rep(1, 
        NR1), rep(2, NR2)))
    data <- as.matrix(counts)
    CD <- new("countData", data = data, replicates = conds, groups = groups)
    CD@libsizes <- getLibsizes(CD)
    CDP.NBML <- getPriors.NB(CD, samplesize = 10000, estimation = "QL", 
        cl = cl)
    CDPost.NBML <- getLikelihoods.NB(CDP.NBML, pET = "BIC", cl = cl)
    all <- topCounts(CDPost.NBML, group = "DE", number = nrow(counts))
    good <- all[all$FDR < cutoff, ]
    
    # Preparing results to return
    result <- new("Result")
    result@data <- all[order(as.numeric(all$FDR)), ]
    result@id <- rownames(result@data)
    result@pval <- result@data$FDR
    stopCluster(cl)
    return(result)
}

run_baySeq_uqn <- function(counts, conds, cutoff, n, runID) {
    
    cl <- makeCluster(2, "SOCK")
    
    # Preparing variables for baySeq run
    NR1 = length(conds[conds == "N"])
    NR2 = length(conds[conds == "T"])
    groups <- list(NDE = c(rep(1, NR1 + NR2)), DE = c(rep(1, 
        NR1), rep(2, NR2)))
    CD <- new("countData", data = as.matrix(UQNnormalization(counts)), 
        replicates = conds, groups = groups)
    
    # CD@libsizes <- getLibsizes(CD);
    meanLib <- mean(apply(CD@data, 2, sum))
    libsizes <- rep(meanLib, times = dim(CD@data)[2])
    names(libsizes) <- colnames(CD@data)
    CD@libsizes <- libsizes    
    CDP.NBML <- getPriors.NB(CD, samplesize = 10000, estimation = "QL", 
        cl = cl)
    CDPost.NBML <- getLikelihoods.NB(CDP.NBML, pET = "BIC", cl = cl)
    all <- topCounts(CDPost.NBML, group = "DE", number = nrow(counts))
    good <- all[all$FDR < cutoff, ]
    
    # Preparing results to return
    result <- new("Result")
    result@data <- all[order(as.numeric(all$FDR)), ]
    result@id <- rownames(result@data)
    result@pval <- result@data$FDR    
    stopCluster(cl)
    return(result)
}


run_baySeq_Mode <- function(counts, conds, cutoff, n, runID, 
    winSize) {
    
    cl <- makeCluster(2, "SOCK")
    # Preparing variables for baySeq run
    NR1 = length(conds[conds == "N"])
    NR2 = length(conds[conds == "T"])
    groups <- list(NDE = c(rep(1, NR1 + NR2)), DE = c(rep(1, 
        NR1), rep(2, NR2)))
    CD <- new("countData", data = as.matrix(normalizeData(counts, 
        conds, runID, winSize)), replicates = conds, groups = groups)
    
    # CD@libsizes <- getLibsizes(CD);
    meanLib <- mean(apply(CD@data, 2, sum))
    libsizes <- rep(meanLib, times = dim(CD@data)[2])
    names(libsizes) <- colnames(CD@data)
    CD@libsizes <- libsizes    
    CDP.NBML <- getPriors.NB(CD, samplesize = 10000, estimation = "QL", 
        cl = cl)
    CDPost.NBML <- getLikelihoods.NB(CDP.NBML, pET = "BIC", cl = cl)
    all <- topCounts(CDPost.NBML, group = "DE", number = nrow(counts))
    good <- all[all$FDR < cutoff, ]
    
    # Preparing results to return
    result <- new("Result")
    result@data <- all[order(as.numeric(all$FDR)), ]
    result@id <- rownames(result@data)
    result@pval <- result@data$FDR
    
    stopCluster(cl)
    return(result)
}


run_baySeq_nde <- function(counts, DElist, conds, cutoff, n, 
    runID) {
    
    cl <- makeCluster(2, "SOCK")
    
    # Preparing variables for baySeq run
    NR1 = length(conds[conds == "N"])
    NR2 = length(conds[conds == "T"])
    groups <- list(NDE = c(rep(1, NR1 + NR2)), DE = c(rep(1, 
        NR1), rep(2, NR2)))
    CD <- new("countData", data = as.matrix(normalizeNDE(counts, 
        DElist, runID)), replicates = conds, groups = groups)
    
    # CD@libsizes <- getLibsizes(CD);
    meanLib <- mean(apply(CD@data, 2, sum))
    libsizes <- rep(meanLib, times = dim(CD@data)[2])
    names(libsizes) <- colnames(CD@data)
    CD@libsizes <- libsizes
    
    CDP.NBML <- getPriors.NB(CD, samplesize = 10000, estimation = "QL", 
        cl = cl)
    CDPost.NBML <- getLikelihoods.NB(CDP.NBML, pET = "BIC", cl = cl)
    all <- topCounts(CDPost.NBML, group = "DE", number = nrow(counts))
    good <- all[all$FDR < cutoff, ]
    
    # Preparing results to return
    result <- new("Result")
    result@data <- all[order(as.numeric(all$FDR)), ]
    result@id <- rownames(result@data)
    result@pval <- result@data$FDR
    
    stopCluster(cl)
    return(result)
}
