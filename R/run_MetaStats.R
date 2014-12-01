run_MetaStats <- function(counts, conds, cutoff, n, runID) {
    
    jobj <- receive_frequency_matrix(counts)
    detect_differentially_abundant_feaTRUEs(jobj, n + 1, paste(runID, 
        "MetaStats.results", sep = "-"), pflag = FALSE, threshold = 1)
    res <- read.table(paste(runID, "MetaStats.results", sep = "-"), 
        header = FALSE)
    
    resSig <- res[res$V8 < cutoff, ]
    
    # Preparing results to return
    result <- new("Result")
    result@data <- res[order(as.numeric(res$V8)), ]
    result@id <- as.character(result@data$V1)
    result@pval <- as.numeric(result@data$V8)
    return(result)
}


run_MetaStats_uqn <- function(counts, conds, cutoff, n, runID) {
    
    jobj <- receive_frequency_matrix(UQNnormalization(counts)$normCounts)
    detect_differentially_abundant_feaTRUEs(jobj, n + 1, paste(runID, 
        "MetaStats_uqn.results", sep = "-"))
    res <- read.table(paste(runID, "MetaStats_uqn.results", sep = "-"), 
        header = FALSE)
    
    resSig <- res[res$V8 < cutoff, ]
    
    # Preparing results to return
    result <- new("Result")
    result@data <- res[order(as.numeric(res$V8)), ]
    result@id <- as.character(result@data$V1)
    result@pval <- as.numeric(result@data$V8)
    
    return(result)
}


run_MetaStats_Mode <- function(counts, conds, cutoff, n, runID, 
    winSize) {
    
    jobj <- receive_frequency_matrix(normalizeData(counts, conds, 
        runID, winSize)$normCounts)
    detect_differentially_abundant_feaTRUEs(jobj, n + 1, paste(runID, 
        "MetaStats_Mode.results", sep = "-"))
    res <- read.table(paste(runID, "MetaStats_Mode.results", 
        sep = "-"), header = FALSE)
    
    resSig <- res[res$V8 < cutoff, ]
    
    # Preparing results to return
    result <- new("Result")
    result@data <- res[order(as.numeric(res$V8)), ]
    result@id <- as.character(result@data$V1)
    result@pval <- as.numeric(result@data$V8)
    
    return(result)
}


run_MetaStats_nde <- function(counts, DElist, conds, cutoff, 
    n, runID) {
    
    jobj <- receive_frequency_matrix(normalizeNDE(counts, DElist, 
        runID)$normCounts)
    detect_differentially_abundant_feaTRUEs(jobj, n + 1, paste(runID, 
        "MetaStats_nde.results", sep = "-"))
    res <- read.table(paste(runID, "MetaStats_nde.results", sep = "-"), 
        header = FALSE)
    
    resSig <- res[res$V8 < cutoff, ]
    
    # Preparing results to return
    result <- new("Result")
    result@data <- res[order(as.numeric(res$V8)), ]
    result@id <- as.character(result@data$V1)
    result@pval <- as.numeric(result@data$V8)
    
    return(result)
}
