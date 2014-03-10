plotROC <- function(obj, DE.methods = c("Cuffdiff", "DESeq", 
    "baySeq", "edgeR", "MetaStats", "NOISeq"), nor.methods = c("default", 
    "Mode", "UQN", "NDE"), plot_type = "o", plot_pch = 20, plot_lwd = 1.75, 
    plot_cex = 1) {
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
    
    plot(NULL, type = plot_type, pch = plot_pch, lwd = plot_lwd, 
        cex = plot_cex, col = "red", xlim = c(0:1), ylim = c(0:1), 
        xlab = "False Positive Rate", ylab = "True Positive Rate")
    legend_text <- c()
    legend_col <- c()
    
    
    
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
    
    # DESeq
    if (is.element("DESeq", method.list) == TRUE) {
        if (class(DESeq)[1] == "Result") {
            print(paste("DESeq AUC: ", DESeq@auc, sep = " "))
            points(DESeq@roc, type = plot_type, pch = plot_pch, 
                lwd = plot_lwd, cex = plot_cex, col = "red4")
            legend_text <- c(legend_text, "DESeq")
            legend_col <- c(legend_col, "red4")
        }
    }
    
    # DESeq_uqn
    if (is.element("DESeq_uqn", method.list) == TRUE) {
        if (class(DESeq_uqn)[1] == "Result") {
            print(paste("DESeq_uqn AUC: ", DESeq_uqn@auc, sep = " "))
            points(DESeq_uqn@roc, type = plot_type, pch = plot_pch, 
                lwd = plot_lwd, cex = plot_cex, col = "red3")
            legend_text <- c(legend_text, "DESeq_uqn")
            legend_col <- c(legend_col, "red3")
        }
    }
    
    # DESeq_Mode
    if (is.element("DESeq_Mode", method.list) == TRUE) {
        if (class(DESeq_Mode)[1] == "Result") {
            print(paste("DESeq_Mode AUC: ", DESeq_Mode@auc, sep = " "))
            points(DESeq_Mode@roc, type = plot_type, pch = plot_pch, 
                lwd = plot_lwd, cex = plot_cex, col = "red2")
            legend_text <- c(legend_text, "DESeq_Mode")
            legend_col <- c(legend_col, "red2")
        }
    }
    
    # DESeq_nde
    if (is.element("DESeq_nde", method.list) == TRUE) {
        if (class(DESeq_nde)[1] == "Result") {
            print(paste("DESeq_nde AUC: ", DESeq_nde@auc, sep = " "))
            points(DESeq_nde@roc, type = plot_type, pch = plot_pch, 
                lwd = plot_lwd, cex = plot_cex, col = "red")
            legend_text <- c(legend_text, "DESeq_nde")
            legend_col <- c(legend_col, "red")
        }
    }
    
    # edgeR
    if (is.element("edgeR", method.list) == TRUE) {
        if (class(edgeR)[1] == "Result") {
            print(paste("edgeR AUC: ", edgeR@auc, sep = " "))
            points(edgeR@roc, type = plot_type, pch = plot_pch, 
                lwd = plot_lwd, cex = plot_cex, col = "green4")
            legend_text <- c(legend_text, "edgeR")
            legend_col <- c(legend_col, "green4")
        }
    }
    
    # edgeR_uqn
    if (is.element("edgeR_uqn", method.list) == TRUE) {
        if (class(edgeR_uqn)[1] == "Result") {
            print(paste("edgeR_uqn AUC: ", edgeR_uqn@auc, sep = " "))
            points(edgeR_uqn@roc, type = plot_type, pch = plot_pch, 
                lwd = plot_lwd, cex = plot_cex, col = "green3")
            legend_text <- c(legend_text, "edgeR_uqn")
            legend_col <- c(legend_col, "green3")
        }
    }
    
    # edgeR_Mode
    if (is.element("edgeR_Mode", method.list) == TRUE) {
        if (class(edgeR_Mode)[1] == "Result") {
            print(paste("edgeR_Mode AUC: ", edgeR_Mode@auc, sep = " "))
            points(edgeR_Mode@roc, type = plot_type, pch = plot_pch, 
                lwd = plot_lwd, cex = plot_cex, col = "green2")
            legend_text <- c(legend_text, "edgeR_Mode")
            legend_col <- c(legend_col, "green2")
        }
    }
    
    # edgeR_nde
    if (is.element("edgeR_nde", method.list) == TRUE) {
        if (class(edgeR_nde)[1] == "Result") {
            print(paste("edgeR_nde AUC: ", edgeR_nde@auc, sep = " "))
            points(edgeR_nde@roc, type = plot_type, pch = plot_pch, 
                lwd = plot_lwd, cex = plot_cex, col = "green")
            legend_text <- c(legend_text, "edgeR_nde")
            legend_col <- c(legend_col, "green")
        }
    }
    
    # baySeq
    if (is.element("baySeq", method.list) == TRUE) {
        if (class(baySeq)[1] == "Result") {
            print(paste("baySeq AUC: ", baySeq@auc, sep = " "))
            points(baySeq@roc, type = plot_type, pch = plot_pch, 
                lwd = plot_lwd, cex = plot_cex, col = "blue")
            legend_text <- c(legend_text, "baySeq")
            legend_col <- c(legend_col, "blue")
        }
    }
    
    # baySeq_uqn
    if (is.element("baySeq_uqn", method.list) == TRUE) {
        if (class(baySeq_uqn)[1] == "Result") {
            print(paste("baySeq_uqn AUC: ", baySeq_uqn@auc, sep = " "))
            points(baySeq_uqn@roc, type = plot_type, pch = plot_pch, 
                lwd = plot_lwd, cex = plot_cex, col = "dodgerblue2")
            legend_text <- c(legend_text, "baySeq_uqn")
            legend_col <- c(legend_col, "dodgerblue2")
        }
    }
    
    # baySeq_Mode
    if (is.element("baySeq_Mode", method.list) == TRUE) {
        if (class(baySeq_Mode)[1] == "Result") {
            print(paste("baySeq_Mode AUC: ", baySeq_Mode@auc, 
                sep = " "))
            points(baySeq_Mode@roc, type = plot_type, pch = plot_pch, 
                lwd = plot_lwd, cex = plot_cex, col = "dodgerblue4")
            legend_text <- c(legend_text, "baySeq_Mode")
            legend_col <- c(legend_col, "dodgerblue4")
        }
    }
    
    # baySeq_nde
    if (is.element("baySeq_nde", method.list) == TRUE) {
        if (class(baySeq_nde)[1] == "Result") {
            print(paste("baySeq_nde AUC: ", baySeq_nde@auc, sep = " "))
            points(baySeq_nde@roc, type = plot_type, pch = plot_pch, 
                lwd = plot_lwd, cex = plot_cex, col = "dodgerblue3")
            legend_text <- c(legend_text, "baySeq_nde")
            legend_col <- c(legend_col, "dodgerblue3")
        }
    }
    
    # NOISeq
    if (is.element("NOISeq", method.list) == TRUE) {
        if (class(NOISeq)[1] == "Result") {
            print(paste("NOISeq AUC: ", NOISeq@auc, sep = " "))
            points(NOISeq@roc, type = plot_type, pch = plot_pch, 
                lwd = plot_lwd, cex = plot_cex, col = "tan4")
            legend_text <- c(legend_text, "NOISeq")
            legend_col <- c(legend_col, "tan4")
        }
    }
    
    # NOISeq_uqn
    if (is.element("NOISeq_uqn", method.list) == TRUE) {
        if (class(NOISeq_uqn)[1] == "Result") {
            print(paste("NOISeq_uqn AUC: ", NOISeq_uqn@auc, sep = " "))
            points(NOISeq_uqn@roc, type = plot_type, pch = plot_pch, 
                lwd = plot_lwd, cex = plot_cex, col = "tomato")
            legend_text <- c(legend_text, "NOISeq_uqn")
            legend_col <- c(legend_col, "tomato")
        }
    }
    
    # NOISeq_Mode
    if (is.element("NOISeq_Mode", method.list) == TRUE) {
        if (class(NOISeq_Mode)[1] == "Result") {
            print(paste("NOISeq_Mode AUC: ", NOISeq_Mode@auc, 
                sep = " "))
            points(NOISeq_Mode@roc, type = plot_type, pch = plot_pch, 
                lwd = plot_lwd, cex = plot_cex, col = "tan1")
            legend_text <- c(legend_text, "NOISeq_Mode")
            legend_col <- c(legend_col, "tan1")
        }
    }
    
    # NOISeq_nde
    if (is.element("NOISeq_nde", method.list) == TRUE) {
        if (class(NOISeq_nde)[1] == "Result") {
            print(paste("NOISeq_nde AUC: ", NOISeq_nde@auc, sep = " "))
            points(NOISeq_nde@roc, type = plot_type, pch = plot_pch, 
                lwd = plot_lwd, cex = plot_cex, col = "tan2")
            legend_text <- c(legend_text, "NOISeq_nde")
            legend_col <- c(legend_col, "tan2")
        }
    }
    
    # Cuffdiff
    if (is.element("Cuffdiff", method.list) == TRUE) {
        if (class(Cuffdiff)[1] == "Result") {
            print(paste("Cuffdiff AUC: ", Cuffdiff@auc, sep = " "))
            points(Cuffdiff@roc, type = plot_type, pch = plot_pch, 
                lwd = plot_lwd, cex = plot_cex, col = "azure4")
            legend_text <- c(legend_text, "Cuffdiff")
            legend_col <- c(legend_col, "azure4")
        }
    }
    
    # Cuffdiff_uqn
    if (is.element("Cuffdiff_uqn", method.list) == TRUE) {
        if (class(Cuffdiff_uqn)[1] == "Result") {
            print(paste("Cuffdiff_uqn AUC: ", Cuffdiff_uqn@auc, 
                sep = " "))
            points(Cuffdiff_uqn@roc, type = plot_type, pch = plot_pch, 
                lwd = plot_lwd, cex = plot_cex, col = "cornflowerblue")
            legend_text <- c(legend_text, "Cuffdiff_uqn")
            legend_col <- c(legend_col, "cornflowerblue")
        }
    }
    
    # Cuffdiff_Mode
    if (is.element("Cuffdiff_Mode", method.list) == TRUE) {
        if (class(Cuffdiff_Mode)[1] == "Result") {
            print(paste("Cuffdiff_Mode AUC: ", Cuffdiff_Mode@auc, 
                sep = " "))
            points(Cuffdiff_Mode@roc, type = plot_type, pch = plot_pch, 
                lwd = plot_lwd, cex = plot_cex, col = "purple")
            legend_text <- c(legend_text, "Cuffdiff_Mode")
            legend_col <- c(legend_col, "purple")
        }
    }
    
    # Cuffdiff_nde
    if (is.element("Cuffdiff_nde", method.list) == TRUE) {
        if (class(Cuffdiff_nde)[1] == "Result") {
            print(paste("Cuffdiff_nde AUC: ", Cuffdiff_nde@auc, 
                sep = " "))
            points(Cuffdiff_nde@roc, type = plot_type, pch = plot_pch, 
                lwd = plot_lwd, cex = plot_cex, col = "purple2")
            legend_text <- c(legend_text, "Cuffdiff_nde")
            legend_col <- c(legend_col, "purple2")
        }
    }
    
    # MetaStats
    if (is.element("MetaStats", method.list) == TRUE) {
        if (class(MetaStats)[1] == "Result") {
            print(paste("MetaStats AUC: ", MetaStats@auc, sep = " "))
            points(MetaStats@roc, type = plot_type, pch = plot_pch, 
                lwd = plot_lwd, cex = plot_cex, col = "lightsteelblue4")
            legend_text <- c(legend_text, "MetaStats")
            legend_col <- c(legend_col, "lightsteelblue4")
        }
    }
    
    # MetaStats_uqn
    if (is.element("MetaStats_uqn", method.list) == TRUE) {
        if (class(MetaStats_uqn)[1] == "Result") {
            print(paste("MetaStats_uqn AUC: ", MetaStats_uqn@auc, 
                sep = " "))
            points(MetaStats_uqn@roc, type = plot_type, pch = plot_pch, 
                lwd = plot_lwd, cex = plot_cex, col = "lightsteelblue2")
            legend_text <- c(legend_text, "MetaStats_uqn")
            legend_col <- c(legend_col, "lightsteelblue2")
        }
    }
    
    # MetaStats_Mode
    if (is.element("MetaStats_Mode", method.list) == TRUE) {
        if (class(MetaStats_Mode)[1] == "Result") {
            print(paste("MetaStats_Mode AUC: ", MetaStats_Mode@auc, 
                sep = " "))
            points(MetaStats_Mode@roc, type = plot_type, pch = plot_pch, 
                lwd = plot_lwd, cex = plot_cex, col = "lightsteelblue")
            legend_text <- c(legend_text, "MetaStats_Mode")
            legend_col <- c(legend_col, "lightsteelblue")
        }
    }
    
    # MetaStats_nde
    if (is.element("MetaStats_nde", method.list) == TRUE) {
        if (class(MetaStats_nde)[1] == "Result") {
            print(paste("MetaStats_nde AUC: ", MetaStats_nde@auc, 
                sep = " "))
            points(MetaStats_nde@roc, type = plot_type, pch = plot_pch, 
                lwd = plot_lwd, cex = plot_cex, col = "lightsteelblue3")
            legend_text <- c(legend_text, "MetaStats_nde")
            legend_col <- c(legend_col, "lightsteelblue3")
        }
    }
    
    legend("bottomright", legend_text, col = legend_col, cex = 1, 
        lty = 1, lwd = 4)
    
}  # end plotROC
