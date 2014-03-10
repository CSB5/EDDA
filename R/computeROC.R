computeROC_matrix <- function(result, cutoff, DElist, NDElist) {
    
    if (is.null(result@data)) {
        result@de <- "NA"
    } else {
        result@de <- nrow(result@data[result@pval < cutoff, ])
    }
    
    ind <- match(DElist, result@id)
    DEE <- rep(0, times = length(result@id))
    DEE[ind] <- 1
    
    pred <- prediction(1 - result@pval, DEE)
    
    perf <- performance(pred, "tpr", "fpr")
    result@roc <- matrix(nrow = length(perf@y.values[[1]]), ncol = 2, 
        data = NA)
    colnames(result@roc) <- c("FPR", "TPR")
    result@roc[, 1] <- perf@x.values[[1]]
    result@roc[, 2] <- perf@y.values[[1]]
    
    perf <- performance(pred, "prec", "rec")
    result@prc <- matrix(nrow = length(perf@y.values[[1]]), ncol = 2, 
        data = NA)
    colnames(result@prc) <- c("Recall", "Precision")
    result@prc[, 1] <- perf@x.values[[1]]
    result@prc[, 2] <- perf@y.values[[1]]
    
    perf <- performance(pred, "auc")
    result@auc <- perf@y.values[[1]]
    return(result)
}  # end computeROC_matrix


computeROC_list <- function(result, cutoff, DElist, NDElist) {
    if (is.null(result@data)) {
        result@de <- "NA"
    } else {
        result@de <- nrow(result@data)
    }
    
    ind <- match(DElist, result@id)
    DEE <- rep(0, times = length(result@id))
    DEE[ind] <- 1
    
    pred <- prediction(result@pval, DEE)
    
    perf <- performance(pred, "tpr", "fpr")
    result@roc <- matrix(nrow = length(perf@y.values[[1]]), ncol = 2, 
        data = NA)
    colnames(result@roc) <- c("FPR", "TPR")
    result@roc[, 1] <- perf@x.values[[1]]
    result@roc[, 2] <- perf@y.values[[1]]
    
    perf <- performance(pred, "prec", "rec")
    result@prc <- matrix(nrow = length(perf@y.values[[1]]), ncol = 2, 
        data = NA)
    colnames(result@prc) <- c("Recall", "Precision")
    result@prc[, 1] <- perf@x.values[[1]]
    result@prc[, 2] <- perf@y.values[[1]]
    
    perf <- performance(pred, "auc")
    result@auc <- perf@y.values[[1]]
    return(result)
}  # end computeROC_list
