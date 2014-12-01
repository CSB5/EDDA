run_DESeq <- function(counts, conds, cutoff, n, runID) {
    # Preparing variables for DESeq run
    cds <- newCountDataSet(counts, conds)
    cds <- estimateSizeFactors(cds)
    if (n > 1) {
        cdsnew <- tryCatch(estimateDispersions(cds, method = "per-condition", 
            sharingMode = "maximum", fitType = "local"), error = function(e) NULL)
        if (is.null(cdsnew)) {
            cdsnew <- tryCatch(estimateDispersions(cds, method = "pooled", 
                sharingMode = "maximum", fitType = "local"), 
                error = function(e) NULL)
        }
        if (is.null(cdsnew)) {
            cdsnew <- tryCatch(estimateDispersions(cds, method = "per-condition", 
                sharingMode = "maximum", fitType = "parametric"), 
                error = function(e) NULL)
        }
        if (is.null(cdsnew)) {
            cds <- estimateDispersions(cds, method = "pooled", 
                sharingMode = "maximum", fitType = "parametric")
        }
        cds <- cdsnew
    } else {
        cdsnew <- tryCatch(estimateDispersions(cds, method = "blind", 
            sharingMode = "fit-only", fitType = "parametric"), 
            error = function(e) NULL)
        if (is.null(cdsnew)) {
            cdsnew <- tryCatch(estimateDispersions(cds, method = "blind", 
                sharingMode = "maximum", fitType = "parametric"), 
                error = function(e) NULL)
        }
        if (is.null(cdsnew)) {
            cdsnew <- tryCatch(estimateDispersions(cds, method = "blind", 
                sharingMode = "fit-only", fitType = "local"), 
                error = function(e) NULL)
        }
        if (is.null(cdsnew)) {
            cds <- estimateDispersions(cds, method = "blind", 
                sharingMode = "maximum", fitType = "local")
        }
        cds <- cdsnew
    }
    res <- nbinomTest(cds, "N", "T")
    res[is.na(res)] <- 0
    resSig <- res[res$pval < cutoff, ]
    
    # Preparing results to return
    result <- new("Result")
    result@data <- res[order(as.numeric(res$pval)), ]
    result@id <- result@data$id
    result@pval <- result@data$pval
    return(result)
}


run_DESeq_uqn <- function(counts, conds, cutoff, n, runID) {
    # Preparing variables for DESeq run
    cds <- newCountDataSet(UQNnormalization(counts)$normCounts, conds)
    sizeFactors(cds) <- c(rep(1, length(conds)))
    if (n > 1) {
        cdsnew <- tryCatch(estimateDispersions(cds, method = "per-condition", 
            sharingMode = "maximum", fitType = "local"), error = function(e) NULL)
        if (is.null(cdsnew)) {
            cdsnew <- tryCatch(estimateDispersions(cds, method = "pooled", 
                sharingMode = "maximum", fitType = "local"), 
                error = function(e) NULL)
        }
        if (is.null(cdsnew)) {
            cdsnew <- tryCatch(estimateDispersions(cds, method = "per-condition", 
                sharingMode = "maximum", fitType = "parametric"), 
                error = function(e) NULL)
        }
        if (is.null(cdsnew)) {
            cds <- estimateDispersions(cds, method = "pooled", 
                sharingMode = "maximum", fitType = "parametric")
        }
        cds <- cdsnew
    } else {
        cdsnew <- tryCatch(estimateDispersions(cds, method = "blind", 
            sharingMode = "fit-only", fitType = "parametric"), 
            error = function(e) NULL)
        if (is.null(cdsnew)) {
            cdsnew <- tryCatch(estimateDispersions(cds, method = "blind", 
                sharingMode = "maximum", fitType = "parametric"), 
                error = function(e) NULL)
        }
        if (is.null(cdsnew)) {
            cdsnew <- tryCatch(estimateDispersions(cds, method = "blind", 
                sharingMode = "fit-only", fitType = "local"), 
                error = function(e) NULL)
        }
        if (is.null(cdsnew)) {
            cds <- estimateDispersions(cds, method = "blind", 
                sharingMode = "maximum", fitType = "local")
        }
        cds <- cdsnew
    }
    res <- nbinomTest(cds, "N", "T")
    res[is.na(res)] <- 0
    resSig <- res[res$pval < cutoff, ]
    
    # Preparing results to return
    result <- new("Result")
    result@data <- res[order(as.numeric(res$pval)), ]
    result@id <- result@data$id
    result@pval <- result@data$pval
    
    return(result)
}


run_DESeq_Mode <- function(counts, conds, cutoff, n, runID, winSize) {
    # Preparing variables for DESeq run
    cds <- newCountDataSet(normalizeData(counts, conds, runID, 
        winSize)$normCounts, conds)
    sizeFactors(cds) <- c(rep(1, length(conds)))
    # cds <- newCountDataSet(counts,conds); sizeFactors(cds) <-
    # computeNormalization(runID, winSizePercentage, minReads);
    if (n > 1) {
        cdsnew <- tryCatch(estimateDispersions(cds, method = "per-condition", 
            sharingMode = "maximum", fitType = "local"), error = function(e) NULL)
        if (is.null(cdsnew)) {
            cdsnew <- tryCatch(estimateDispersions(cds, method = "pooled", 
                sharingMode = "maximum", fitType = "local"), 
                error = function(e) NULL)
        }
        if (is.null(cdsnew)) {
            cdsnew <- tryCatch(estimateDispersions(cds, method = "per-condition", 
                sharingMode = "maximum", fitType = "parametric"), 
                error = function(e) NULL)
        }
        if (is.null(cdsnew)) {
            cds <- estimateDispersions(cds, method = "pooled", 
                sharingMode = "maximum", fitType = "parametric")
        }
        cds <- cdsnew
    } else {
        cdsnew <- tryCatch(estimateDispersions(cds, method = "blind", 
            sharingMode = "fit-only", fitType = "parametric"), 
            error = function(e) NULL)
        if (is.null(cdsnew)) {
            cdsnew <- tryCatch(estimateDispersions(cds, method = "blind", 
                sharingMode = "maximum", fitType = "parametric"), 
                error = function(e) NULL)
        }
        if (is.null(cdsnew)) {
            cdsnew <- tryCatch(estimateDispersions(cds, method = "blind", 
                sharingMode = "fit-only", fitType = "local"), 
                error = function(e) NULL)
        }
        if (is.null(cdsnew)) {
            cds <- estimateDispersions(cds, method = "blind", 
                sharingMode = "maximum", fitType = "local")
        }
        cds <- cdsnew
    }
    res <- nbinomTest(cds, "N", "T")
    res[is.na(res)] <- 0
    resSig <- res[res$pval < cutoff, ]
    
    # Preparing results to return
    result <- new("Result")
    result@data <- res[order(as.numeric(res$pval)), ]
    result@id <- result@data$id
    result@pval <- result@data$pval
    
    return(result)
}


run_DESeq_nde <- function(counts, DElist, conds, cutoff, n, runID) {
    # Preparing variables for DESeq run
    cds <- newCountDataSet(normalizeNDE(counts, DElist, runID)$normCounts, 
        conds)
    sizeFactors(cds) <- c(rep(1, n * 2))
    # cds <- newCountDataSet(counts,conds); sizeFactors(cds) <-
    # computeNormalization(runID, winSizePercentage, minReads);
    if (n > 1) {
        cdsnew <- tryCatch(estimateDispersions(cds, method = "per-condition", 
            sharingMode = "maximum", fitType = "local"), error = function(e) NULL)
        if (is.null(cdsnew)) {
            cdsnew <- tryCatch(estimateDispersions(cds, method = "pooled", 
                sharingMode = "maximum", fitType = "local"), 
                error = function(e) NULL)
        }
        if (is.null(cdsnew)) {
            cdsnew <- tryCatch(estimateDispersions(cds, method = "per-condition", 
                sharingMode = "maximum", fitType = "parametric"), 
                error = function(e) NULL)
        }
        if (is.null(cdsnew)) {
            cds <- estimateDispersions(cds, method = "pooled", 
                sharingMode = "maximum", fitType = "parametric")
        }
        cds <- cdsnew
    } else {
        cdsnew <- tryCatch(estimateDispersions(cds, method = "blind", 
            sharingMode = "fit-only", fitType = "parametric"), 
            error = function(e) NULL)
        if (is.null(cdsnew)) {
            cdsnew <- tryCatch(estimateDispersions(cds, method = "blind", 
                sharingMode = "maximum", fitType = "parametric"), 
                error = function(e) NULL)
        }
        if (is.null(cdsnew)) {
            cdsnew <- tryCatch(estimateDispersions(cds, method = "blind", 
                sharingMode = "fit-only", fitType = "local"), 
                error = function(e) NULL)
        }
        if (is.null(cdsnew)) {
            cds <- estimateDispersions(cds, method = "blind", 
                sharingMode = "maximum", fitType = "local")
        }
        cds <- cdsnew
    }
    res <- nbinomTest(cds, "N", "T")
    res[is.na(res)] <- 0
    resSig <- res[res$pval < cutoff, ]
    
    # Preparing results to return
    result <- new("Result")
    result@data <- res[order(as.numeric(res$pval)), ]
    result@id <- result@data$id
    result@pval <- result@data$pval
    
    
    return(result)
}
