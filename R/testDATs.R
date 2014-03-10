testDATs <- function(data, numCores = 10, minCountsThreshold = 0, 
    DE.methods = c("Cuffdiff", "DESeq", "baySeq", "edgeR", "MetaStats", 
        "NOISeq"), nor.methods = c("default", "Mode", "UQN", 
        "NDE"), method.list = NULL) {
    cutoff = 1
    winSize = 10
    if (.Platform$OS.type == "windows") 
        numCores = 1
    result <- list(data = data, DESeq = NULL, DESeq_uqn = NULL, 
        DESeq_Mode = NULL, DESeq_nde = NULL, edgeR = NULL, edgeR_uqn = NULL, 
        edgeR_Mode = NULL, edgeR_nde = NULL, baySeq = NULL, baySeq_uqn = NULL, 
        baySeq_Mode = NULL, baySeq_nde = NULL, NOISeq = NULL, 
        NOISeq_uqn = NULL, NOISeq_Mode = NULL, NOISeq_nde = NULL, 
        Cuffdiff = NULL, Cuffdiff_uqn = NULL, Cuffdiff_Mode = NULL, 
        Cuffdiff_nde = NULL, MetaStats = NULL, MetaStats_uqn = NULL, 
        MetaStats_Mode = NULL, MetaStats_nde = NULL, filterCounts = NULL)
    if (is.null(method.list)) {
        method.list = NULL
        method.num = 1
        for (i in 1:length(DE.methods)) {
            for (j in 1:length(nor.methods)) {
                if (nor.methods[j] == "Mode") 
                  nor.methods[j] = "Mode"
                if (nor.methods[j] == "default") {
                  method.list[method.num] = DE.methods[i]
                } else {
                  method.list[method.num] = paste(DE.methods[i], 
                    nor.methods[j], sep = "_")
                }
                method.num = method.num + 1
            }
        }
    }
    # Step 1: Read countsFile
    counts <- filterGenesByCounts(data.frame(data$count), minCountsThreshold)
    result$filterCounts <- counts
    NR1 = length(data$dataLabel[data$dataLabel == 0])
    NR2 = length(data$dataLabel[data$dataLabel == 1])
    n <- NR1
    runID = ""
    # Setting up variables for DESeq test
    conds <- c(rep("N", NR1), rep("T", NR2))
    analyses <- list(job_DESeq = function() if (is.element("DESeq", 
        method.list) == TRUE) {
        print("Starting DESeq analysis...")
        DESeq <- tryCatch(run_DESeq(counts, conds, cutoff, n, 
            runID), error = function(e) e)
        if (class(DESeq)[1] != "Result") {
            print("Error running DESeq:")
            print(DESeq$message)
        }
        print("DESeq analysis completed.")
        return(DESeq)
    }, job_DESeq_uqn = function() if (is.element("DESeq_UQN", 
        method.list) == TRUE) {
        print("Starting DESeq_uqn analysis...")
        DESeq_uqn <- tryCatch(run_DESeq_uqn(counts, conds, cutoff, 
            n, runID), error = function(e) e)
        if (class(DESeq_uqn)[1] != "Result") {
            print("Error running DESeq_uqn:")
            print(DESeq_uqn$message)
        }
        print("DESeq_uqn analysis completed.")
        return(DESeq_uqn)
    }, job_DESeq_Mode = function() if (is.element("DESeq_Mode", 
        method.list) == TRUE) {
        print("Starting DESeq_Mode analysis...")
        DESeq_Mode <- tryCatch(run_DESeq_Mode(counts, conds, 
            cutoff, n, runID, winSize), error = function(e) e)
        if (class(DESeq_Mode)[1] != "Result") {
            print("Error running DESeq_Mode:")
            print(DESeq_Mode$message)
        }
        print("DESeq_Mode analysis completed.")
        return(DESeq_Mode)
    }, job_DESeq_nde = function() if (is.element("DESeq_NDE", 
        method.list) == TRUE) {
        print("Starting DESeq_nde analysis...")
        DESeq_nde <- tryCatch(run_DESeq_nde(counts, data$DiffAbundList$geneName, 
            conds, cutoff, n, runID), error = function(e) e)
        if (class(DESeq_nde)[1] != "Result") {
            print("Error running DESeq_nde:")
            print(DESeq_nde$message)
        }
        print("DESeq_nde analysis completed.")
        return(DESeq_nde)
    }, job_edgeR = function() if (is.element("edgeR", method.list) == 
        TRUE) {
        print("Starting edgeR analysis...")
        edgeR <- tryCatch(run_edgeR(counts, conds, cutoff, n, 
            runID), error = function(e) e)
        if (class(edgeR)[1] != "Result") {
            print("Error running edgeR:")
            print(edgeR$message)
        }
        print("edgeR analysis completed.")
        return(edgeR)
    }, job_edgeR_uqn = function() if (is.element("edgeR_UQN", 
        method.list) == TRUE) {
        print("Starting edgeR_uqn analysis...")
        edgeR_uqn <- tryCatch(run_edgeR_uqn(counts, conds, cutoff, 
            n, runID), error = function(e) e)
        if (class(edgeR_uqn)[1] != "Result") {
            print("Error running edgeR_uqn:")
            print(edgeR_uqn$message)
        }
        print("edgeR_uqn analysis completed.")
        return(edgeR_uqn)
    }, job_edgeR_Mode = function() if (is.element("edgeR_Mode", 
        method.list) == TRUE) {
        print("Starting edgeR_Mode analysis...")
        edgeR_Mode <- tryCatch(run_edgeR_Mode(counts, conds, 
            cutoff, n, runID, winSize), error = function(e) e)
        if (class(edgeR_Mode)[1] != "Result") {
            print("Error running edgeR_Mode:")
            print(edgeR_Mode$message)
        }
        print("edgeR_Mode analysis completed.")
        return(edgeR_Mode)
    }, job_edgeR_nde = function() if (is.element("edgeR_NDE", 
        method.list) == TRUE) {
        print("Starting edgeR_nde analysis...")
        edgeR_nde <- tryCatch(run_edgeR_nde(counts, data$DiffAbundList$geneName, 
            conds, cutoff, n, runID), error = function(e) e)
        if (class(edgeR_nde)[1] != "Result") {
            print("Error running edgeR_nde:")
            print(edgeR_nde$message)
        }
        print("edgeR_nde analysis completed.")
        return(edgeR_nde)
    }, job_baySeq = function() if (is.element("baySeq", method.list) == 
        TRUE) {
        print("Starting baySeq analysis...")
        baySeq <- tryCatch(run_baySeq(counts, conds, cutoff, 
            n, runID), error = function(e) e)
        if (class(baySeq)[1] != "Result") {
            print("Error running baySeq:")
            print(baySeq$message)
        }
        print("baySeq analysis completed.")
        return(baySeq)
    }, job_NOISeq = function() if (is.element("NOISeq", method.list) == 
        TRUE) {
        print("Starting NOISeq analysis...")
        NOISeq <- tryCatch(run_NOISeq(counts, conds, cutoff, 
            n, runID), error = function(e) e)
        if (class(NOISeq)[1] != "Result") {
            print("Error running NOISeq:")
            print(NOISeq$message)
        }
        print("NOISeq analysis completed.")
        return(NOISeq)
    }, job_NOISeq_uqn = function() if (is.element("NOISeq_UQN", 
        method.list) == TRUE) {
        print("Starting NOISeq_uqn analysis...")
        NOISeq_uqn <- tryCatch(run_NOISeq_uqn(counts, conds, 
            cutoff, n, runID), error = function(e) e)
        if (class(NOISeq_uqn)[1] != "Result") {
            print("Error running NOISeq_uqn:")
            print(NOISeq_uqn$message)
        }
        print("NOISeq_uqn analysis completed.")
        return(NOISeq_uqn)
    }, job_NOISeq_Mode = function() if (is.element("NOISeq_Mode", 
        method.list) == TRUE) {
        print("Starting NOISeq_Mode analysis...")
        NOISeq_Mode <- tryCatch(run_NOISeq_Mode(counts, conds, 
            cutoff, n, runID, winSize), error = function(e) e)
        if (class(NOISeq_Mode)[1] != "Result") {
            print("Error running NOISeq_Mode:")
            print(NOISeq_Mode$message)
        }
        print("NOISeq_Mode analysis completed.")
        return(NOISeq_Mode)
    }, job_NOISeq_nde = function() if (is.element("NOISeq_NDE", 
        method.list) == TRUE) {
        print("Starting NOISeq_nde analysis...")
        NOISeq_nde <- tryCatch(run_NOISeq_nde(counts, 
            data$DiffAbundList$geneName, 
            conds, cutoff, n, runID), error = function(e) e)
        if (class(NOISeq_nde)[1] != "Result") {
            print("Error running NOISeq_nde:")
            print(NOISeq_nde$message)
        }
        print("NOISeq_nde analysis completed.")
        return(NOISeq_nde)
    }, job_Cuffdiff = function() if (is.element("Cuffdiff", method.list) == 
        TRUE) {
        print("Starting Cuffdiff analysis...")
        Cuffdiff <- tryCatch(run_Cuffdiff(counts, conds, cutoff, 
            n, runID), error = function(e) e)
        if (class(Cuffdiff)[1] != "Result") {
            print("Error running Cuffdiff:")
            print(Cuffdiff$message)
        }
        print("Cuffdiff analysis completed.")
        return(Cuffdiff)
    }, job_Cuffdiff_uqn = function() if (is.element("Cuffdiff_UQN", 
        method.list) == TRUE) {
        print("Starting Cuffdiff_uqn analysis...")
        Cuffdiff_uqn <- tryCatch(run_Cuffdiff_uqn(counts, conds, 
            cutoff, n, runID), error = function(e) e)
        if (class(Cuffdiff_uqn)[1] != "Result") {
            print("Error running Cuffdiff_uqn:")
            print(Cuffdiff_uqn$message)
        }
        print("Cuffdiff_uqn analysis completed.")
        return(Cuffdiff_uqn)
    }, job_Cuffdiff_Mode = function() if (is.element("Cuffdiff_Mode", 
        method.list) == TRUE) {
        print("Starting Cuffdiff_Mode analysis...")
        Cuffdiff_Mode <- tryCatch(run_Cuffdiff_Mode(counts, conds, 
            cutoff, n, runID, winSize), error = function(e) e)
        if (class(Cuffdiff_Mode)[1] != "Result") {
            print("Error running Cuffdiff_Mode:")
            print(Cuffdiff_Mode$message)
        }
        print("Cuffdiff_Mode analysis completed.")
        return(Cuffdiff_Mode)
    }, job_Cuffdiff_nde = function() if (is.element("Cuffdiff_NDE", 
        method.list) == TRUE) {
        print("Starting Cuffdiff_nde analysis...")
        Cuffdiff_nde <- tryCatch(run_Cuffdiff_nde(counts, 
            data$DiffAbundList$geneName, 
            conds, cutoff, n, runID), error = function(e) e)
        if (class(Cuffdiff_nde)[1] != "Result") {
            print("Error running Cuffdiff_nde:")
            print(Cuffdiff_nde$message)
        }
        print("Cuffdiff_nde analysis completed.")
        return(Cuffdiff_nde)
    }, job_MetaStats = function() if (is.element("MetaStats", 
        method.list) == TRUE) {
        print("Starting MetaStats analysis...")
        MetaStats <- tryCatch(run_MetaStats(counts, conds, cutoff, 
            n, runID), error = function(e) e)
        if (class(MetaStats)[1] != "Result") {
            print("Error running MetaStats:")
            print(MetaStats$message)
        }
        print("MetaStats analysis completed.")
        return(MetaStats)
    }, job_MetaStats_uqn = function() if (is.element("MetaStats_UQN", 
        method.list) == TRUE) {
        print("Starting MetaStats_uqn analysis...")
        MetaStats_uqn <- tryCatch(run_MetaStats_uqn(counts, conds, 
            cutoff, n, runID), error = function(e) e)
        if (class(MetaStats_uqn)[1] != "Result") {
            print("Error running MetaStats_uqn:")
            print(MetaStats_uqn$message)
        }
        print("MetaStats_uqn analysis completed.")
        return(MetaStats_uqn)
    }, job_MetaStats_Mode = function() if (is.element("MetaStats_Mode", 
        method.list) == TRUE) {
        print("Starting MetaStats_Mode analysis...")
        MetaStats_Mode <- tryCatch(run_MetaStats_Mode(counts, 
            conds, cutoff, n, runID, winSize), error = function(e) e)
        if (class(MetaStats_Mode)[1] != "Result") {
            print("Error running MetaStats_Mode:")
            print(MetaStats_Mode$message)
        }
        print("MetaStats_Mode analysis completed.")
        return(MetaStats_Mode)
    }, job_MetaStats_nde = function() if (is.element("MetaStats_NDE", 
        method.list) == TRUE) {
        print("Starting MetaStats_nde analysis...")
        MetaStats_nde <- tryCatch(run_MetaStats_nde(counts, 
            data$DiffAbundList$geneName, 
            conds, cutoff, n, runID), error = function(e) e)
        if (class(MetaStats_nde)[1] != "Result") {
            print("Error running MetaStats_nde:")
            print(MetaStats_nde$message)
        }
        print("MetaStats_nde analysis completed.")
        return(MetaStats_nde)
    })
    parallel <- mclapply(analyses, function(f) f(), mc.cores = numCores)
    result$DESeq <- parallel$job_DESeq
    result$DESeq_uqn <- parallel$job_DESeq_uqn
    result$DESeq_Mode <- parallel$job_DESeq_Mode
    result$DESeq_nde <- parallel$job_DESeq_nde
    result$edgeR <- parallel$job_edgeR
    result$edgeR_uqn <- parallel$job_edgeR_uqn
    result$edgeR_Mode <- parallel$job_edgeR_Mode
    result$edgeR_nde <- parallel$job_edgeR_nde
    result$baySeq <- parallel$job_baySeq
    result$NOISeq <- parallel$job_NOISeq
    result$NOISeq_uqn <- parallel$job_NOISeq_uqn
    result$NOISeq_Mode <- parallel$job_NOISeq_Mode
    result$NOISeq_nde <- parallel$job_NOISeq_nde
    result$Cuffdiff <- parallel$job_Cuffdiff
    result$Cuffdiff_uqn <- parallel$job_Cuffdiff_uqn
    result$Cuffdiff_Mode <- parallel$job_Cuffdiff_Mode
    result$Cuffdiff_nde <- parallel$job_Cuffdiff_nde
    result$MetaStats <- parallel$job_MetaStats
    result$MetaStats_uqn <- parallel$job_MetaStats_uqn
    result$MetaStats_Mode <- parallel$job_MetaStats_Mode
    result$MetaStats_nde <- parallel$job_MetaStats_nde
    # Step 2.10: Running baySeq_uqn analysis
    if (is.element("baySeq_UQN", method.list) == TRUE) {
        print("Starting baySeq_uqn analysis...")
        baySeq_uqn <- tryCatch(run_baySeq_uqn(counts, conds, 
            cutoff, n, runID), error = function(e) e)
        if (class(baySeq_uqn)[1] != "Result") {
            print("Error running baySeq_uqn:")
            print(baySeq_uqn$message)
        }
        print("baySeq_uqn analysis completed.")
        result$baySeq_uqn <- baySeq_uqn
    }
    # Step 2.11: Running baySeq_Mode analysis
    if (is.element("baySeq_Mode", method.list) == TRUE) {
        print("Starting baySeq_Mode analysis...")
        baySeq_Mode <- tryCatch(run_baySeq_Mode(counts, conds, 
            cutoff, n, runID, winSize), error = function(e) e)
        if (class(baySeq_Mode)[1] != "Result") {
            print("Error running baySeq_Mode:")
            print(baySeq_Mode$message)
        }
        print("baySeq_Mode analysis completed.")
        result$baySeq_Mode <- baySeq_Mode
    }
    # Step 2.12: Running baySeq_nde analysis
    if (is.element("baySeq_NDE", method.list) == TRUE) {
        print("Starting baySeq_nde analysis...")
        baySeq_nde <- tryCatch(run_baySeq_nde(counts, 
            data$DiffAbundList$geneName, 
            conds, cutoff, n, runID), error = function(e) e)
        if (class(baySeq_nde)[1] != "Result") {
            print("Error running baySeq_nde:")
            print(baySeq_nde$message)
        }
        print("baySeq_nde analysis completed.")
        result$baySeq_nde <- baySeq_nde
    }
    return(result)
}
