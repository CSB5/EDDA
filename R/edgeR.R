# Overloading methods from edgeR

exactTestMode <- function(object, pair = 1:2, dispersion = "auto", 
    rejection.region = "doubletail", big.count = 900, prior.count.total = 0.5) 
    # Calculates exact p-values for the differential expression
# levels of tags in the two groups being compared.  Davis
# McCarthy, Gordon Smyth.  Created September 2009. Last
# modified 8 July 2012.
{
    # Check input
    if (!is(object, "DGEList")) 
        stop("Currently only supports DGEList objects \n\tas the object argument.")
    if (length(pair) != 2) 
        stop("Pair must be of length 2.")
    rejection.region <- match.arg(rejection.region, c("doubletail", 
        "deviance", "smallp"))
    
    # Get group names
    group <- as.factor(object$samples$group)
    levs.group <- levels(group)
    if (is.numeric(pair)) 
        pair <- levs.group[pair] else pair <- as.character(pair)
    if (!all(pair %in% levs.group)) 
        stop("At least one element of given pair \n\tis not a group.\n Groups are: ", 
            paste(levs.group, collapse = " "))
    
    # Get dispersion vector
    if (is.null(dispersion)) 
        dispersion <- "auto"
    if (is.character(dispersion)) {
        dispersion <- match.arg(dispersion, c("auto", "common", 
            "trended", "tagwise"))
        dispersion <- switch(dispersion, common = object$common.dispersion, 
            trended = object$trended.dispersion, 
            tagwise = object$tagwise.dispersion, 
            auto = getDispersion(object))
        if (is.null(dispersion)) 
            stop("specified dispersion not found in object")
        if (is.na(dispersion[1])) 
            stop("dispersion is NA")
    }
    ldisp <- length(dispersion)
    ntags <- nrow(object$counts)
    if (ldisp != 1 && ldisp != ntags) 
        stop("Dispersion provided by user must have \n\tlength either 1 or 
        	the number of tags in the DGEList object.")
    if (ldisp == 1) 
        dispersion <- rep(dispersion, ntags)
    
    # Reduce to two groups
    group <- as.character(group)
    j <- group %in% pair
    y <- object$counts[, j, drop = FALSE]
    lib.size <- object$samples$lib.size[j]
    norm.factors <- object$samples$norm.factors[j]
    group <- group[j]
    if (is.null(rownames(y))) 
        rownames(y) <- paste("tag", 1:ntags, sep = ".")
    
    # Normalized library sizes
    lib.size <- lib.size * norm.factors
    offset <- log(lib.size)
    lib.size.average <- exp(mean(offset))
    
    # Average abundance
    abundance <- mglmOneGroup(y, dispersion = dispersion, offset = offset)
    logCPM <- (abundance + log(1e+06))/log(2)
    
    # logFC
    prior.count <- lib.size
    prior.count <- prior.count.total * prior.count/sum(prior.count)
    j1 <- group == pair[1]
    n1 <- sum(j1)
    if (n1 == 0) 
        stop("No libraries for", pair[1])
    y1 <- y[, j1, drop = FALSE]
    abundance1 <- mglmOneGroup(y1 + matrix(prior.count[j1], ntags, 
        n1, byrow = TRUE), offset = offset[j1])
    j2 <- group == pair[2]
    n2 <- sum(j2)
    if (n1 == 0) 
        stop("No libraries for", pair[2])
    y2 <- y[, j2, drop = FALSE]
    abundance2 <- mglmOneGroup(y2 + matrix(prior.count[j2], ntags, 
        n2, byrow = TRUE), offset = offset[j2])
    logFC <- (abundance2 - abundance1)/log(2)
    
    # Equalize library sizes
    e <- exp(abundance)
    input.mean <- matrix(e, ntags, n1)
    output.mean <- input.mean * lib.size.average
    input.mean <- t(t(input.mean) * lib.size[j1])
    input.mean <- matrix(e, ntags, n2)
    output.mean <- input.mean * lib.size.average
    input.mean <- t(t(input.mean) * lib.size[j2])
    exact.pvals <- switch(rejection.region, doubletail = exactTestDoubleTail(y1, 
        y2, dispersion = dispersion, big.count = big.count), 
        deviance = exactTestByDeviance(y1, y2, dispersion = dispersion, 
            big.count = big.count), smallp = exactTestBySmallP(y1, 
            y2, dispersion = dispersion, big.count = big.count))
    
    de.out <- data.frame(logFC = logFC, logCPM = logCPM, PValue = exact.pvals)
    rn <- rownames(object$counts)
    if (!is.null(rn)) 
        rownames(de.out) <- make.unique(rn)
    new("DGEExact", list(table = de.out, comparison = pair, genes = object$genes))
}

equalizeLibSizesMode <- function(object, dispersion = 0, common.lib.size = NULL) 
# Normalize counts so that the library sizes can be treated
# as equal.  Uses a quantile-to-quantile transformation so
# that new count counts are equivalent deviates on the
# equalized scale.  Davis McCarthy, Gordon Smyth.  Created
# July 2009. Last modified 25 July 2012.
{
    d <- dim(object)
    ntags <- d[1]
    nlibs <- d[2]
    lib.size <- object$samples$lib.size * object$samples$norm.factors
    if (is.null(common.lib.size)) 
        common.lib.size <- exp(mean(log(lib.size)))
    levs.group <- unique(object$samples$group)
    if (length(dispersion) == 1) 
        dispersion <- rep(dispersion, ntags)
    
    input.mean <- output.mean <- matrix(0, ntags, nlibs)
    for (i in 1:length(levs.group)) {
        j <- object$samples$group == levs.group[i]
        beta <- mglmOneGroup(object$counts[, j, drop = FALSE], 
            dispersion = dispersion, offset = log(lib.size[j]))
        lambda <- exp(beta)
        input.mean[, j] <- matrix(lambda, ncol = 1) %*% matrix(lib.size[j], 
            nrow = 1)
        output.mean[, j] <- matrix(lambda, ncol = 1) %*% matrix(common.lib.size, 
            nrow = 1, ncol = sum(j))
    }
    # pseudo <- q2qnbinom(object$counts, input.mean=input.mean,
    # output.mean=output.mean,dispersion=disp)
    # pseudo[pseudo<0]<-0
    list(pseudo.counts = object$counts, common.lib.size = common.lib.size)
}

estimateCommonDispMode <- function(object, tol = 1e-06, rowsum.filter = 5, 
    verbose = FALSE) 
# Do two iterations of calculating pseudodata and estimating
# common dispersion Davis McCarthy, Mark Robinson, Gordon
# Smyth.  Created 2009. Last modified 2 Aug 2012.
{
    if (!is(object, "DGEList")) 
        stop("Currently supports DGEList objects")
    group <- object$samples$group <- as.factor(object$samples$group)
    
    if (all(tabulate(group) <= 1)) {
        warning("There is no replication, setting dispersion to NA.")
        object$common.dispersion <- NA
        return(object)
    }
    
    tags.used <- rowSums(object$counts) > rowsum.filter
    pseudo.obj <- object[tags.used, ]
    
    # Start from small dispersion
    disp <- 0.01
    for (i in 1:2) {
        out <- equalizeLibSizesMode(object, dispersion = disp)
        pseudo.obj$counts <- out$pseudo.counts[tags.used, , drop = FALSE]
        y <- splitIntoGroups(pseudo.obj)
        delta <- optimize(commonCondLogLikDerDelta, interval = c(1e-04, 
            100/(100 + 1)), tol = tol, maximum = TRUE, y = y, 
            der = 0)
        delta <- delta$maximum
        disp <- delta/(1 - delta)
    }
    if (verbose) 
        cat("Disp =", round(disp, 5), ", BCV =", round(sqrt(disp), 
            4), "\n")
    object$common.dispersion <- disp
    object$pseudo.counts <- out$pseudo.counts
    
    # Average logCPM
    effective.lib.size <- object$samples$lib.size * object$samples$norm.factors
    abundance <- mglmOneGroup(object$counts, dispersion = disp, 
        offset = log(effective.lib.size))
    object$logCPM <- log1p(exp(abundance + log(1e+06)))/log(2)
    object$pseudo.lib.size <- out$common.lib.size
    object
}
