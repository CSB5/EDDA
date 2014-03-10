computeAUC <- function(obj, cutoff = 1, numCores = 10, 
    DE.methods = c("Cuffdiff", "DESeq", "baySeq", "edgeR", "MetaStats",
     "NOISeq"), nor.methods = c("default", "Mode", "UQN", "NDE")) {
    if (.Platform$OS.type == "windows") 
        numCores = 1
    DElist = obj$data$DiffAbundList$geneName
    
    DElist <- intersect(DElist, rownames(obj$filterCounts))
    
    NDElist = setdiff(rownames(obj$filterCounts), DElist)
    
    method.list = NULL
    method.num = 1
    
    for (i in 1:length(DE.methods)) {
        for (j in 1:length(nor.methods)) {
            if (nor.methods[j] == "Mode") 
                nor.methods[j] = "Mode"
            if (nor.methods[j] == "UQN") 
                nor.methods[j] = "uqn"
            if (nor.methods[j] == "NDE") 
                nor.methods[j] = "nde"
            if (nor.methods[j] == "default") {
                method.list[method.num] = DE.methods[i]
            } else {
                method.list[method.num] = paste(DE.methods[i], 
                  nor.methods[j], sep = "_")
            }
            method.num = method.num + 1
        }
    }
    
    DESeq <- obj[["DESeq"]]
    DESeq_uqn <- obj$DESeq_uqn
    DESeq_Mode <- obj$DESeq_Mode
    DESeq_nde <- obj$DESeq_nde
    edgeR <- obj[["edgeR"]]
    edgeR_uqn <- obj$edgeR_uqn
    edgeR_Mode <- obj$edgeR_Mode
    edgeR_nde <- obj$edgeR_nde
    baySeq <- obj[["baySeq"]]
    baySeq_uqn <- obj$baySeq_uqn
    baySeq_Mode <- obj$baySeq_Mode
    baySeq_nde <- obj$baySeq_nde
    NOISeq <- obj[["NOISeq"]]
    NOISeq_uqn <- obj$NOISeq_uqn
    NOISeq_Mode <- obj$NOISeq_Mode
    NOISeq_nde <- obj$NOISeq_nde
    Cuffdiff <- obj[["Cuffdiff"]]
    Cuffdiff_uqn <- obj$Cuffdiff_uqn
    Cuffdiff_Mode <- obj$Cuffdiff_Mode
    Cuffdiff_nde <- obj$Cuffdiff_nde
    MetaStats <- obj[["MetaStats"]]
    MetaStats_uqn <- obj$MetaStats_uqn
    MetaStats_Mode <- obj$MetaStats_Mode
    MetaStats_nde <- obj$MetaStats_nde
    
    plots <- list(job_DESeq = function() if (is.element("DESeq", 
        method.list) == TRUE) {
        if (class(DESeq)[1] == "Result") {
            print("Computing plot values for DESeq...")
            DESeq <- computeROC_matrix(DESeq, cutoff, DElist, 
                NDElist)
            print("Computation for DESeq completed.")
            return(DESeq)
        }
    }, job_DESeq_uqn = function() if (is.element("DESeq_uqn", 
        method.list) == TRUE) {
        if (class(DESeq_uqn)[1] == "Result") {
            print("Computing plot values for DESeq_uqn...")
            DESeq_uqn <- computeROC_matrix(DESeq_uqn, cutoff, 
                DElist, NDElist)
            print("Computation for DESeq_uqn completed.")
            return(DESeq_uqn)
        }
    }, job_DESeq_Mode = function() if (is.element("DESeq_Mode", 
        method.list) == TRUE) {
        if (class(DESeq_Mode)[1] == "Result") {
            print("Computing plot values for DESeq_Mode...")
            DESeq_Mode <- computeROC_matrix(DESeq_Mode, cutoff, 
                DElist, NDElist)
            print("Computation for DESeq_Mode completed.")
            return(DESeq_Mode)
        }
    }, job_DESeq_nde = function() if (is.element("DESeq_nde", 
        method.list) == TRUE) {
        if (class(DESeq_nde)[1] == "Result") {
            print("Computing plot values for DESeq_nde...")
            DESeq_nde <- computeROC_matrix(DESeq_nde, cutoff, 
                DElist, NDElist)
            print("Computation for DESeq_nde completed.")
            return(DESeq_nde)
        }
    }, job_edgeR = function() if (is.element("edgeR", method.list) == 
        TRUE) {
        if (class(edgeR)[1] == "Result") {
            print("Computing plot values for edgeR...")
            edgeR <- computeROC_matrix(edgeR, cutoff, DElist, 
                NDElist)
            print("Computation for edgeR completed.")
            return(edgeR)
        }
    }, job_edgeR_uqn = function() if (is.element("edgeR_uqn", 
        method.list) == TRUE) {
        if (class(edgeR_uqn)[1] == "Result") {
            print("Computing plot values for edgeR_uqn...")
            edgeR_uqn <- computeROC_matrix(edgeR_uqn, cutoff, 
                DElist, NDElist)
            print("Computation for edgeR_uqn completed.")
            return(edgeR_uqn)
        }
    }, job_edgeR_Mode = function() if (is.element("edgeR_Mode", 
        method.list) == TRUE) {
        if (class(edgeR_Mode)[1] == "Result") {
            print("Computing plot values for edgeR_Mode...")
            edgeR_Mode <- computeROC_matrix(edgeR_Mode, cutoff, 
                DElist, NDElist)
            print("Computation for edgeR_Mode completed.")
            return(edgeR_Mode)
        }
    }, job_edgeR_nde = function() if (is.element("edgeR_nde", 
        method.list) == TRUE) {
        if (class(edgeR_nde)[1] == "Result") {
            print("Computing plot values for edgeR_nde...")
            edgeR_nde <- computeROC_matrix(edgeR_nde, cutoff, 
                DElist, NDElist)
            print("Computation for edgeR_nde completed.")
            return(edgeR_nde)
        }
    }, job_baySeq = function() if (is.element("baySeq", method.list) == 
        TRUE) {
        if (class(baySeq)[1] == "Result") {
            print("Computing plot values for baySeq...")
            baySeq <- computeROC_matrix(baySeq, cutoff, DElist, 
                NDElist)
            print("Computation for baySeq completed.")
            return(baySeq)
        }
    }, job_baySeq_uqn = function() if (is.element("baySeq_uqn", 
        method.list) == TRUE) {
        if (class(baySeq_uqn)[1] == "Result") {
            print("Computing plot values for baySeq_uqn...")
            baySeq_uqn <- computeROC_matrix(baySeq_uqn, cutoff, 
                DElist, NDElist)
            print("Computation for baySeq_uqn completed.")
            return(baySeq_uqn)
        }
    }, job_baySeq_Mode = function() if (is.element("baySeq_Mode", 
        method.list) == TRUE) {
        if (class(baySeq_Mode)[1] == "Result") {
            print("Computing plot values for baySeq_Mode...")
            baySeq_Mode <- computeROC_matrix(baySeq_Mode, cutoff, 
                DElist, NDElist)
            print("Computation for baySeq_Mode completed.")
            return(baySeq_Mode)
        }
    }, job_baySeq_nde = function() if (is.element("baySeq_nde", 
        method.list) == TRUE) {
        if (class(baySeq_nde)[1] == "Result") {
            print("Computing plot values for baySeq_nde...")
            baySeq_nde <- computeROC_matrix(baySeq_nde, cutoff, 
                DElist, NDElist)
            print("Computation for baySeq_nde completed.")
            return(baySeq_nde)
        }
    }, job_NOISeq = function() if (is.element("NOISeq", method.list) == 
        TRUE) {
        if (class(NOISeq)[1] == "Result") {
            print("Computing plot values for NOISeq...")
            ## source(NOISeq_src);
            NOISeq <- computeROC_list(NOISeq, cutoff, DElist, 
                NDElist)
            print("Computation for NOISeq completed.")
            return(NOISeq)
        }
    }, job_NOISeq_uqn = function() if (is.element("NOISeq_uqn", 
        method.list) == TRUE) {
        if (class(NOISeq_uqn)[1] == "Result") {
            print("Computing plot values for NOISeq_uqn...")
            ## source(NOISeq_src);
            NOISeq_uqn <- computeROC_list(NOISeq_uqn, cutoff, 
                DElist, NDElist)
            print("Computation for NOISeq_uqn completed.")
            return(NOISeq_uqn)
        }
    }, job_NOISeq_Mode = function() if (is.element("NOISeq_Mode", 
        method.list) == TRUE) {
        if (class(NOISeq_Mode)[1] == "Result") {
            print("Computing plot values for NOISeq_Mode...")
            ## source(NOISeq_src);
            NOISeq_Mode <- computeROC_list(NOISeq_Mode, cutoff, 
                DElist, NDElist)
            print("Computation for NOISeq_Mode completed.")
            return(NOISeq_Mode)
        }
    }, job_NOISeq_nde = function() if (is.element("NOISeq_nde", 
        method.list) == TRUE) {
        if (class(NOISeq_nde)[1] == "Result") {
            print("Computing plot values for NOISeq_nde...")
            ## source(NOISeq_src);
            NOISeq_nde <- computeROC_list(NOISeq_nde, cutoff, 
                DElist, NDElist)
            print("Computation for NOISeq_nde completed.")
            return(NOISeq_nde)
        }
    }, job_Cuffdiff = function() if (is.element("Cuffdiff", method.list) == 
        TRUE) {
        if (class(Cuffdiff)[1] == "Result") {
            print("Computing plot values for Cuffdiff...")
            Cuffdiff <- computeROC_matrix(Cuffdiff, cutoff, DElist, 
                NDElist)
            print("Computation for Cuffdiff completed.")
            return(Cuffdiff)
        }
    }, job_Cuffdiff_uqn = function() if (is.element("Cuffdiff_uqn", 
        method.list) == TRUE) {
        if (class(Cuffdiff_uqn)[1] == "Result") {
            print("Computing plot values for Cuffdiff_uqn...")
            Cuffdiff_uqn <- computeROC_matrix(Cuffdiff_uqn, cutoff, 
                DElist, NDElist)
            print("Computation for Cuffdiff_uqn completed.")
            return(Cuffdiff_uqn)
        }
    }, job_Cuffdiff_Mode = function() if (is.element("Cuffdiff_Mode", 
        method.list) == TRUE) {
        if (class(Cuffdiff_Mode)[1] == "Result") {
            print("Computing plot values for Cuffdiff_Mode...")
            Cuffdiff_Mode <- computeROC_matrix(Cuffdiff_Mode, 
                cutoff, DElist, NDElist)
            print("Computation for Cuffdiff_Mode completed.")
            return(Cuffdiff_Mode)
        }
    }, job_Cuffdiff_nde = function() if (is.element("Cuffdiff_nde", 
        method.list) == TRUE) {
        if (class(Cuffdiff_nde)[1] == "Result") {
            print("Computing plot values for Cuffdiff_nde...")
            Cuffdiff_nde <- computeROC_matrix(Cuffdiff_nde, cutoff, 
                DElist, NDElist)
            print("Computation for Cuffdiff_nde completed.")
            return(Cuffdiff_nde)
        }
    }, job_MetaStats = function() if (is.element("MetaStats", 
        method.list) == TRUE) {
        if (class(MetaStats)[1] == "Result") {
            print("Computing plot values for MetaStats...")
            MetaStats <- computeROC_matrix(MetaStats, cutoff, 
                DElist, NDElist)
            print("Computation for MetaStats completed.")
            return(MetaStats)
        }
    }, job_MetaStats_uqn = function() if (is.element("MetaStats_uqn", 
        method.list) == TRUE) {
        if (class(MetaStats_uqn)[1] == "Result") {
            print("Computing plot values for MetaStats_uqn...")
            MetaStats_uqn <- computeROC_matrix(MetaStats_uqn, 
                cutoff, DElist, NDElist)
            print("Computation for MetaStats_uqn completed.")
            return(MetaStats_uqn)
        }
    }, job_MetaStats_Mode = function() if (is.element("MetaStats_Mode", 
        method.list) == TRUE) {
        if (class(MetaStats_Mode)[1] == "Result") {
            print("Computing plot values for MetaStats_Mode...")
            MetaStats_Mode <- computeROC_matrix(MetaStats_Mode, 
                cutoff, DElist, NDElist)
            print("Computation for MetaStats_Mode completed.")
            return(MetaStats_Mode)
        }
    }, job_MetaStats_nde = function() if (is.element("MetaStats_nde", 
        method.list) == TRUE) {
        if (class(MetaStats_nde)[1] == "Result") {
            print("Computing plot values for MetaStats_nde...")
            MetaStats_nde <- computeROC_matrix(MetaStats_nde, 
                cutoff, DElist, NDElist)
            print("Computation for MetaStats_nde completed.")
            return(MetaStats_nde)
        }
    })
    
    parallel <- mclapply(plots, function(f) f(), mc.cores = numCores)
    
    obj$DESeq <- parallel$job_DESeq
    obj$DESeq_uqn <- parallel$job_DESeq_uqn
    obj$DESeq_Mode <- parallel$job_DESeq_Mode
    obj$DESeq_nde <- parallel$job_DESeq_nde
    obj$edgeR <- parallel$job_edgeR
    obj$edgeR_uqn <- parallel$job_edgeR_uqn
    obj$edgeR_Mode <- parallel$job_edgeR_Mode
    obj$edgeR_nde <- parallel$job_edgeR_nde
    obj$baySeq <- parallel$job_baySeq
    obj$baySeq_uqn <- parallel$job_baySeq_uqn
    obj$baySeq_Mode <- parallel$job_baySeq_Mode
    obj$baySeq_nde <- parallel$job_baySeq_nde
    obj$NOISeq <- parallel$job_NOISeq
    obj$NOISeq_uqn <- parallel$job_NOISeq_uqn
    obj$NOISeq_Mode <- parallel$job_NOISeq_Mode
    obj$NOISeq_nde <- parallel$job_NOISeq_nde
    obj$Cuffdiff <- parallel$job_Cuffdiff
    obj$Cuffdiff_uqn <- parallel$job_Cuffdiff_uqn
    obj$Cuffdiff_Mode <- parallel$job_Cuffdiff_Mode
    obj$Cuffdiff_nde <- parallel$job_Cuffdiff_nde
    obj$MetaStats <- parallel$job_MetaStats
    obj$MetaStats_uqn <- parallel$job_MetaStats_uqn
    obj$MetaStats_Mode <- parallel$job_MetaStats_Mode
    obj$MetaStats_nde <- parallel$job_MetaStats_nde
    
    return(obj)
    
}
