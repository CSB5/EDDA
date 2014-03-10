run_Cuffdiff <- function(counts, conds, cutoff, n, runID) {
    res <- Cuffdiff(counts, length(conds[conds == "N"]), 
      length(conds[conds == "T"]))
    resSig <- res[res$pValue < cutoff, ]
    
    # Preparing results to return
    result <- new("Result")
    result@data <- res[order(as.numeric(res$pValue)), ]
    result@id <- result@data$id
    result@pval <- as.numeric(result@data$pValue)
    return(result)
}

run_Cuffdiff_uqn <- function(counts, conds, cutoff, n, runID) {
    res <- Cuffdiff_uqn(counts, length(conds[conds == "N"]), 
        length(conds[conds == "T"]))
    resSig <- res[res$pValue < cutoff, ]
    
    # Preparing results to return
    result <- new("Result")
    result@data <- res[order(as.numeric(res$pValue)), ]
    result@id <- result@data$id
    result@pval <- as.numeric(result@data$pValue)
    return(result)
}

run_Cuffdiff_Mode <- function(counts, conds, cutoff, n, runID, 
    winSize) {
    res <- Cuffdiff_Mode(counts, conds, length(conds[conds == 
        "N"]), length(conds[conds == "T"]), runID, winSize)
    resSig <- res[res$pValue < cutoff, ]
    
    # Preparing results to return
    result <- new("Result")
    result@data <- res[order(as.numeric(res$pValue)), ]
    result@id <- result@data$id
    result@pval <- as.numeric(result@data$pValue)
    
    return(result)
}

run_Cuffdiff_nde <- function(counts, DElist, conds, cutoff, n, 
    runID) {
    res <- Cuffdiff_nde(counts, DElist, length(conds[conds == 
        "N"]), length(conds[conds == "T"]), runID)
    resSig <- res[res$pValue < cutoff, ]
    # Preparing results to return
    result <- new("Result")
    result@data <- res[order(as.numeric(res$pValue)), ]
    result@id <- result@data$id
    result@pval <- as.numeric(result@data$pValue)
    
    return(result)
}
