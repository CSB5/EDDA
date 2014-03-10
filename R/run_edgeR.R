run_edgeR <- function(counts, conds, cutoff, n, runID) {
    # Preparing variables for edgeR run
    d <- DGEList(counts = counts, group = conds)
    d = calcNormFactors(d, method = "TMM")
    d = estimateCommonDisp(d)
    de.com = exactTest(d)
    options(digits = 4)
    detags.com = rownames(topTags(de.com)$table)
    all <- topTags(de.com, n = nrow(counts))$table
    good = sum(all$PValue < cutoff)
    goodList = topTags(de.com, n = good)
    
    # Preparing results to return
    result <- new("Result")
    result@data <- all[order(as.numeric(all$PValue)), ]
    result@id <- rownames(result@data)
    result@pval <- result@data$PValue
    result@de <- d$common.dispersion
    
    return(result)
}


run_edgeR_uqn <- function(counts, conds, cutoff, n, runID) {
    d <- DGEList(counts = counts, group = conds)
    d = calcNormFactors(d, method = "upperquartile")
    d = estimateCommonDisp(d)
    de.com = exactTest(d)
    options(digits = 4)
    detags.com = rownames(topTags(de.com)$table)
    all <- topTags(de.com, n = nrow(counts))$table
    good = sum(all$PValue < cutoff)
    goodList = topTags(de.com, n = good)
    
    # Preparing results to return
    result <- new("Result")
    result@data <- all[order(as.numeric(all$PValue)), ]
    result@id <- rownames(result@data)
    result@pval <- result@data$PValue
    
    return(result)
}


run_edgeR_Mode <- function(counts, conds, cutoff, n, runID, winSize) {
    d <- DGEList(counts = normalizeData(counts, conds, runID, 
        winSize), group = conds)
    d$samples$norm.factors <- rep(1, length(conds))
    
    d = estimateCommonDispMode(d)
    de.com = exactTestMode(d)
    options(digits = 4)
    detags.com = rownames(topTags(de.com)$table)
    all <- topTags(de.com, n = nrow(counts))$table
    good = sum(all$PValue < cutoff)
    goodList = topTags(de.com, n = good)
    
    # Preparing results to return
    result <- new("Result")
    result@data <- all[order(as.numeric(all$PValue)), ]
    result@id <- rownames(result@data)
    result@pval <- result@data$PValue
    
    return(result)
}


run_edgeR_nde <- function(counts, DElist, conds, cutoff, n, runID) {
    # Preparing variables for edgeR run
    d <- DGEList(counts = normalizeNDE(counts, DElist, runID), 
        group = conds)
    d$samples$norm.factors <- rep(1, n * 2)
    # d$samples$norm.factors <- normFactors; d =
    # calcNormFactors(d);
    d = estimateCommonDispMode(d)
    de.com = exactTestMode(d)
    options(digits = 4)
    detags.com = rownames(topTags(de.com)$table)
    all <- topTags(de.com, n = nrow(counts))$table
    good = sum(all$PValue < cutoff)
    goodList = topTags(de.com, n = good)
    
    # Preparing results to return
    result <- new("Result")
    result@data <- all[order(as.numeric(all$PValue)), ]
    result@id <- rownames(result@data)
    result@pval <- result@data$PValue
    return(result)
}
