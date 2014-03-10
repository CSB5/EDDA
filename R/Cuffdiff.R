Cuffdiff <- function(counts, cond1, cond2) {
    rpm <- sweep(counts, 2, apply(counts, 2, sum)/1e+06, "/")
    # normalize counts by per million read to get RPM
    results <- data.frame(matrix(nrow = length(rpm[, 1]), ncol = 7))
    # create empty data frame to hold results
    rownames(results) <- rownames(rpm)
    # setting rownames of results matrix
    colnames(results) <- c("id", "RPM_A", "RPMvar_A", "RPM_B", 
        "RPMvar_B", "pValue", "pAdj")
    # setting colnames of results matrix
    
    for (i in seq(from = 1, to = length(rpm[, 1]), by = 1)) {
        id <- rownames(rpm)[i]
        if (cond1 == 1) {
            RPM_A <- rpm[i, 1:cond1]
        } else {
            RPM_A <- rowMeans(rpm[i, 1:cond1])
        }
        RPMvar_A <- RPM_A
        if (cond2 == 1) {
            RPM_B <- rpm[i, (cond1 + 1):(cond1 + cond2)]
        } else {
            RPM_B <- rowMeans(rpm[i, (cond1 + 1):(cond1 + cond2)])
        }
        RPMvar_B <- RPM_B
        pVal <- .Call("cuffdiff_wrapper", RPM_A, RPMvar_A, RPM_B, 
            RPMvar_B)
        results[i, 1:6] <- c(id, RPM_A, RPMvar_A, RPM_B, RPMvar_B, 
            pVal)
        # storing results in table
    }  # end for loop
    results[, 7] <- p.adjust(results[, 6], method = "BH")
    
    return(results)
}


Cuffdiff_uqn <- function(counts, cond1, cond2) {
    uqn <- UQNnormalization(counts)
    results <- data.frame(matrix(nrow = length(uqn[, 1]), ncol = 7))
    # create empty data frame to hold results
    rownames(results) <- rownames(uqn)
    # setting rownames of results matrix
    colnames(results) <- c("id", "UQN_A", "UQNvar_A", "UQN_B", 
        "UQNvar_B", "pValue", "pAdj")
    # setting colnames of results matrix
    
    for (i in seq(from = 1, to = length(uqn[, 1]), by = 1)) {
        id <- rownames(uqn)[i]
        if (cond1 == 1) {
            UQN_A <- uqn[i, 1:cond1]
        } else {
            UQN_A <- mean(uqn[i, 1:cond1])
        }
        UQNvar_A <- UQN_A
        if (cond2 == 1) {
            UQN_B <- uqn[i, (cond1 + 1):(cond1 + cond2)]
        } else {
            UQN_B <- mean(uqn[i, (cond1 + 1):(cond1 + cond2)])
        }
        UQNvar_B <- UQN_B
        pVal <- .Call("cuffdiff_wrapper", UQN_A, UQNvar_A, UQN_B, 
            UQNvar_B)
        results[i, 1:6] <- c(id, UQN_A, UQNvar_A, UQN_B, UQNvar_B, 
            pVal)
        # storing results in table
    }  # end for loop
    results[, 7] <- p.adjust(results[, 6], method = "BH")
    
    return(results)
}


Cuffdiff_Mode <- function(counts, conds, cond1, cond2, runID, 
    winSize) {
    normCounts <- normalizeData(counts, conds, runID, winSize)
    results <- data.frame(matrix(nrow = length(normCounts[, 1]), 
        ncol = 7))
    # create empty data frame to hold results
    rownames(results) <- rownames(normCounts)
    # setting rownames of results matrix
    colnames(results) <- c("id", "counts_A", "var_A", "counts_B", 
        "var_B", "pValue", "pAdj")
    # setting colnames of results matrix
    
    for (i in seq(from = 1, to = length(normCounts[, 1]), by = 1)) {
        id <- rownames(normCounts)[i]
        if (cond1 == 1) {
            counts_A <- normCounts[i, 1:cond1]
        } else {
            counts_A <- mean(normCounts[i, 1:cond1])
        }
        var_A <- counts_A
        if (cond2 == 1) {
            counts_B <- normCounts[i, (cond1 + 1):(cond1 + cond2)]
        } else {
            counts_B <- mean(normCounts[i, (cond1 + 1):(cond1 + 
                cond2)])
        }
        var_B <- counts_B
        pVal <- .Call("cuffdiff_wrapper", counts_A, var_A, counts_B, 
            var_B)
        results[i, 1:6] <- c(id, counts_A, var_A, counts_B, var_B, 
            pVal)
        # storing results in table
    }  # end for loop
    results[, 7] <- p.adjust(results[, 6], method = "BH")
    
    return(results)
}


Cuffdiff_nde <- function(counts, DElist, cond1, cond2, runID) {
    normCounts <- normalizeNDE(counts, DElist, runID)
    results <- data.frame(matrix(nrow = length(normCounts[, 1]), 
        ncol = 7))
    # create empty data frame to hold results
    rownames(results) <- rownames(normCounts)
    # setting rownames of results matrix
    colnames(results) <- c("id", "counts_A", "var_A", "counts_B", 
        "var_B", "pValue", "pAdj")
    # setting colnames of results matrix
    
    for (i in seq(from = 1, to = length(normCounts[, 1]), by = 1)) {
        id <- rownames(normCounts)[i]
        if (cond1 == 1) {
            counts_A <- normCounts[i, 1:cond1]
        } else {
            counts_A <- mean(normCounts[i, 1:cond1])
        }
        var_A <- counts_A
        if (cond2 == 1) {
            counts_B <- normCounts[i, (cond1 + 1):(cond1 + cond2)]
        } else {
            counts_B <- mean(normCounts[i, (cond1 + 1):(cond1 + 
                cond2)])
        }
        var_B <- counts_B
        pVal <- .Call("cuffdiff_wrapper", counts_A, var_A, counts_B, 
            var_B)
        results[i, 1:6] <- c(id, counts_A, var_A, counts_B, var_B, 
            pVal)
        # storing results in table
    }  # end for loop
    results[, 7] <- p.adjust(results[, 6], method = "BH")
    
    return(results)
}
