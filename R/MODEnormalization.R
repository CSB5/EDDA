
############################ Supporting functions ##
setClass(Class = "Result", representation(roc = "matrix", prc = "matrix", 
    de = "numeric", auc = "numeric", data = "data.frame", id = "character", 
    pval = "numeric"))  # end setClass


filterGenesByCounts <- function(counts, minCountsThreshold) {
    rowSum <- apply(counts, 1, sum)
    return(counts[rowSum>minCountsThreshold & !is.na(rowSum), ])
}


MODEnormalization <- function(counts, conds, runID, winSize) {
    # version 4 of Mode normalization. New algorithm which finds
    # three peaks or tightest peak as the NDE set of genes.
    
    oldcounts <- counts
    ######################################### 
    meanExp <- apply(counts, 1, mean)
    counts <- counts[meanExp > 10, ]
    
    ####################################### 
    lib1 <- apply(as.matrix(counts[, conds == "N"]), 2, mean)
    lib2 <- apply(as.matrix(counts[, conds == "T"]), 2, mean)
    
    normfactor <- lib1/median(lib1)
    counts[, conds == "N"] <- round(sweep(as.matrix(counts[, 
        conds == "N"]), 2, normfactor, "/"))
    normfactor <- lib2/median(lib2)
    counts[, conds == "T"] <- round(sweep(as.matrix(counts[, 
        conds == "T"]), 2, normfactor, "/"))
    ####################################### 
    
    
    
    a <- apply(as.matrix(counts[, conds == "N"]), 1, mean)
    b <- apply(as.matrix(counts[, conds == "T"]), 1, mean)
    
    ind <- (a != 0) & (b != 0)
    a <- a[ind]
    b <- b[ind]
    sizefactor <- log2(a/b)
    names(sizefactor) <- names(a)
    sizefactor <- sizefactor[!is.na(sizefactor)]
    
    best_bw <- 0.28
    best_peak <- 3
    
    for (bw in seq(0.5, 0.1, by = -0.02)) {
        dens <- density(sizefactor, bw = bw)
        pk <- dens$x[peaks(dens$y)]
        pky <- dens$y[peaks(dens$y)]
        
        ind <- pky > 0.02
        pk <- pk[ind]
        pky <- pky[ind]
        
        if (length(pk) == 3) {
            if ((pk[2] - pk[1]) >= 0.5 && (pk[3] - pk[2]) >= 
                0.5) {
                best_peak = 3
            } else if ((pk[2] - pk[1]) >= 0.5 || (pk[3] - pk[2]) >= 
                0.5) {
                best_peak = 2
            } else if ((pk[2] - pk[1]) < 0.5 || (pk[3] - pk[2]) < 
                0.5) {
                best_peak = 1
            }
            break
        }
    }
    
    for (bw in seq(0.5, 0.1, by = -0.02)) {
        dens <- density(sizefactor, bw = bw)
        pk <- dens$x[peaks(dens$y)]
        pky <- dens$y[peaks(dens$y)]
        
        ind <- pky > 0.02
        pk <- pk[ind]
        pky <- pky[ind]
        
        if (length(pk) == best_peak) {
            best_bw = bw
            break
        }
    }
    
    dens <- density(sizefactor, bw = best_bw)
    pk <- dens$x[peaks(dens$y)]
    pky <- dens$y[peaks(dens$y)]
    
    ind <- dens$y[peaks(dens$y)] > 0.05
    pk <- pk[ind]
    pky <- pky[ind]
    
    if (length(pk) > 1) {
        pkyod <- order(pky, decreasing = TRUE)
        if (pky[pkyod[1]]/pky[pkyod[2]] > 5) {
            pk = pk[pkyod[1]]
            pky = pky[pkyod[1]]
        }
    }
    
    if (length(pk) == 1) {
        NDEpk = pk[1]
        NDEpky = pky[1]
    } else if (length(pk) == 3) {
        NDEpk = pk[2]
        NDEpky = pky[2]
    } else if (length(pk) == 2) {
        tempy <- dens$y[dens$x > pk[1] & dens$x < pk[2]]
        pky_bottom <- min(tempy)
        
        if (abs(pky[1] - pky_bottom) < 0.05) {
            NDEpk <- pk[2]
            NDEpky <- pky[2]
        } else if (abs(pky[2] - pky_bottom) < 0.05) {
            NDEpk <- pk[1]
            NDEpky <- pky[1]
        } else {
            peakness <- NULL
            lefty <- dens$y[which.min(abs(dens$x - pk[1]))]/2
            tempx <- dens$x[dens$x < pk[1]]
            tempy <- dens$y[dens$x < pk[1]]
            
            peakness <- c(peakness, pk[1] - tempx[which.min(abs(tempy - 
                lefty))])
            
            righty <- dens$y[which.min(abs(dens$x - pk[2]))]/2
            tempx <- dens$x[dens$x > pk[2]]
            tempy <- dens$y[dens$x > pk[2]]
            
            peakness <- c(peakness, tempx[which.min(abs(tempy - 
                righty))] - pk[2])
            
            NDEpk <- pk[which.min(peakness)]
            NDEpky <- pky[which.min(peakness)]
            print("Two modes observed!!!")
        }
    } else {
        pkyod <- order(pky, decreasing = TRUE)[1:3]
        pk = pk[pkyod]
        NDEpk = pk[2]
        NDEpky = pky[pkyod][2]
    }
    
    winSize = NDEpky * 0.1 * dim(counts)[1]
    gg <- names(sort(abs(sizefactor - NDEpk)))[1:winSize]
    
    ggCounts <- counts[gg, ]
    i <- apply(ggCounts == 0, 1, any)
    if (any(i)) 
        ggCounts <- ggCounts[!i, , drop = FALSE]
    
    gg <- names(sort(apply(ggCounts, 1, mean), 
        decreasing = TRUE))[1:max(dim(ggCounts)[1]/2, 1)]
    ggCounts <- oldcounts[gg, ]
    tmp <- exp(colMeans(log(ggCounts)))
    normfactor <- tmp/mean(tmp)
    
    return(normfactor)
}  # end MODEnormalization


normalizeData <- function(counts, conds, runID, winSize) {
    # returns the normalized counts
    normFactors <- MODEnormalization(counts, conds, runID, winSize)
    normCounts <- round(as.matrix(sweep(counts, 2, normFactors, 
        "/")))
    
    return(normCounts)
}  # end normalizeData


normalizeNDE <- function(data, DElist, runID) {
    # normalize by GM of all NDE gene counts
    counts <- data + 1
    # DEfile <- read.delim(paste(runID,'-DEgeneList.txt',
    # sep=''), header=TRUE); DElist <- DEfile$geneName;
    NDElist <- setdiff(rownames(counts), DElist)
    
    tmp <- exp(colMeans(log(counts[NDElist, ])))
    normfactor <- tmp/mean(tmp)
    
    normCounts <- round(as.matrix(sweep(counts, 2, normfactor, 
        "/")))
    
    return(normCounts)
}  # end normalizeData



peaks <- function(series, span = 11) {
    z <- embed(series, span)
    s <- span%/%2
    # v<- max.col(z) == 1 + s
    ind <- apply(z, 1, which.max)
    v <- ind == (1 + s)
    result <- c(rep(FALSE, s), v)
    result <- result[1:(length(result) - s)]
    result
}


UQNnormalization <- function(counts, p = 0.75) {
    data <- counts
    i <- apply(data <= 0, 1, all)
    if (any(i)) 
        data <- data[!i, , drop = FALSE]
    f <- apply(data, 2, function(x) quantile(x, p = p))
    normFactors <- f/exp(mean(log(f)))
    
    normCounts <- round(as.matrix(sweep(data, 2, normFactors, 
        "/")))
    return(normCounts)
}


sinceros <- function(datos, k) {
    # Replacing counts=0 with counts=k
    datos0 <- datos
    if (is.null(k)) {
        mini0 <- min(datos[noceros(datos, num = FALSE, k = 0)])
        kc <- mini0/2
        datos0[datos0 == 0] <- kc
    } else {
        datos0[datos0 == 0] <- k
    }
    datos0
}  # end sinceros

