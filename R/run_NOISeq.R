run_NOISeq <- function(counts, conds, cutoff, n, runID) {    
    # Preparing variables for NOISeq run
    NR1 = length(conds[conds == "N"])
    NR2 = length(conds[conds == "T"])
    mydata <- readData(counts, cond1 = c(1:NR1), cond2 = c((NR1 + 
        1):(NR1 + NR2)))
    rownames(mydata[[1]]) <- rownames(counts)
    rownames(mydata[[2]]) <- rownames(counts)
    
    # Running NOISeq
    if (n > 1) {
        myresults <- noiseq(mydata[[2]], mydata[[1]], repl = "tech", 
            norm = "rpkm", q = 0.8)
        all <- noiseq(mydata[[2]], mydata[[1]], repl = "tech", 
            norm = "rpkm", q = 0)
    } else {
        myresults <- noiseq(mydata[[2]], mydata[[1]], k = NULL, 
            norm = "rpkm", q = 0.8, pnr = 0.2, nss = 5, v = 0.02)
        all <- noiseq(mydata[[2]], mydata[[1]], k = NULL, norm = "rpkm", 
            q = 0, pnr = 0.2, nss = 5, v = 0.02)
    }
    
    # Preparing results to return
    result <- new("Result")
    result@data <- 
      as.data.frame(na.omit(all$probab[order(as.numeric(-all$probab))]))
    result@id <- rownames(result@data)
    result@pval <- as.numeric(result@data[, 1])
    
    return(result)
}  # end run_NOISeq


run_NOISeq_uqn <- function(counts, conds, cutoff, n, runID) {
    
    # Preparing variables for NOISeq run
    NR1 = length(conds[conds == "N"])
    NR2 = length(conds[conds == "T"])
    
    counts <- UQNnormalization(counts)$normCounts
    ind <- apply(counts, 1, mean) > 0
    counts <- counts[ind, ]
    
    mydata <- readData(counts, cond1 = c(1:NR1), cond2 = c((NR1 + 
        1):(NR1 + NR2)))
    rownames(mydata[[1]]) <- rownames(counts)
    rownames(mydata[[2]]) <- rownames(counts)
    
    # Running NOISeq
    if (n > 1) {
        myresults <- noiseq(mydata[[2]], mydata[[1]], repl = "tech", 
            norm = "n", q = 0.8)
        all <- noiseq(mydata[[2]], mydata[[1]], repl = "tech", 
            norm = "n", q = 0)
    } else {
        myresults <- noiseq(mydata[[2]], mydata[[1]], k = NULL, 
            norm = "n", q = 0.8, pnr = 0.2, nss = 5, v = 0.02)
        all <- noiseq(mydata[[2]], mydata[[1]], k = NULL, norm = "n", 
            q = 0, pnr = 0.2, nss = 5, v = 0.02)
    }
    
    # Preparing results to return
    result <- new("Result")
    result@data <- 
      as.data.frame(na.omit(all$probab[order(as.numeric(-all$probab))]))
    result@id <- rownames(result@data)
    result@pval <- as.numeric(result@data[, 1])
    
    return(result)
}  # end run_NOISeq_uqn


run_NOISeq_Mode <- function(counts, conds, cutoff, n, runID, 
    winSize) {
    
    
    # Preparing variables for NOISeq run
    NR1 = length(conds[conds == "N"])
    NR2 = length(conds[conds == "T"])
    
    counts <- normalizeData(counts, conds, runID, winSize)$normCounts
    ind <- apply(counts, 1, mean) > 0
    counts <- counts[ind, ]
    
    mydata <- readData(counts, cond1 = c(1:NR1), cond2 = c((NR1 + 
        1):(NR1 + NR2)))
    rownames(mydata[[1]]) <- rownames(counts)
    rownames(mydata[[2]]) <- rownames(counts)
    
    # Running NOISeq
    if (n > 1) {
        myresults <- noiseq(mydata[[2]], mydata[[1]], repl = "tech", 
            norm = "n", q = 0.8)
        all <- noiseq(mydata[[2]], mydata[[1]], repl = "tech", 
            norm = "n", q = 0)
    } else {
        myresults <- noiseq(mydata[[2]], mydata[[1]], k = NULL, 
            norm = "n", q = 0.8, pnr = 0.2, nss = 5, v = 0.02)
        all <- noiseq(mydata[[2]], mydata[[1]], k = NULL, norm = "n", 
            q = 0, pnr = 0.2, nss = 5, v = 0.02)
    }
    
    # Preparing results to return
    result <- new("Result")
    result@data <- 
      as.data.frame(na.omit(all$probab[order(as.numeric(-all$probab))]))
    result@id <- rownames(result@data)
    result@pval <- as.numeric(result@data[, 1])
    
    
    return(result)
}  # end run_NOISeq_Mode


run_NOISeq_nde <- function(counts, DElist, conds, cutoff, n, 
    runID) {
    
    # Preparing variables for NOISeq run
    NR1 = length(conds[conds == "N"])
    NR2 = length(conds[conds == "T"])
    
    mydata <- readData(normalizeNDE(counts, DElist, runID)$normCounts, cond1 = c(1:NR1), 
        cond2 = c((NR1 + 1):(NR1 + NR2)))
    rownames(mydata[[1]]) <- rownames(counts)
    rownames(mydata[[2]]) <- rownames(counts)
    
    # Running NOISeq
    if (n > 1) {
        myresults <- noiseq(mydata[[2]], mydata[[1]], repl = "tech", 
            norm = "n", q = 0.8)
        all <- noiseq(mydata[[2]], mydata[[1]], repl = "tech", 
            norm = "n", q = 0)
    } else {
        myresults <- noiseq(mydata[[2]], mydata[[1]], k = NULL, 
            norm = "n", q = 0.8, pnr = 0.2, nss = 5, v = 0.02)
        all <- noiseq(mydata[[2]], mydata[[1]], k = NULL, norm = "n", 
            q = 0, pnr = 0.2, nss = 5, v = 0.02)
    }
    
    # Preparing results to return
    result <- new("Result")
    result@data <- 
      as.data.frame(na.omit(all$probab[order(as.numeric(-all$probab))]))
    result@id <- rownames(result@data)
    result@pval <- as.numeric(result@data[, 1])
    
    return(result)
}  # end run_NOISeq_nde


