generateData <- function(SimulModel = "Full", SampleVar = "medium", 
    ControlRep = 5, CaseRep = ControlRep, EntityCount = 1000, 
    FC = "Norm(2,1)", perDiffAbund = 0.1, upPDA = perDiffAbund/2, 
    downPDA = perDiffAbund/2, numDataPoints = 100, AbundProfile = "HBR", modelFile = NULL, 
    minAbund = 10, varLibsizes = 0.1, outlier=FALSE,perOutlier=0.15, factorOutlier=100, inputCount = NULL, inputLabel = NULL, 
    SimulType = "auto") {
    
    #Check setting correct first to avoid meaningless waiting if wrong
    if (SimulModel == "ModelFree" || SimulModel == "ModelFreeMn") {
      if (!is.null(modelFile))
      {
        if (modelFile %in% c("SingleCell")){
          
        }else if (file.exists(modelFile)){
          model.data <- read.delim(modelFile, header = TRUE, stringsAsFactors = FALSE, 
                                   row.names = 1)
          if (dim(model.data)[2]<=1)
          {
            stop("ERROR: No replicates!")
          }
        }else{
          stop("ERROR: modelFile not correct! Set modelFile to NULL or correct name or correct path!")
        }
        
      }else if (is.null(modelFile) && (is.null(inputCount) || is.null(inputLabel))){
        stop("ERROR: using ModelFree approach without model file and pilot data for sampling.")
      }
    }
    
    if (outlier == TRUE && (perOutlier>1 || perOutlier<0))
    {
      stop("ERROR: Outlier parameter should be in range 0 to 1!")
    }
    
    model = SimulModel
    NR1 = ControlRep
    NR2 = CaseRep
    
    EC = EntityCount
    perDE = perDiffAbund
    upDE = upPDA
    downDE = downPDA
    ND = numDataPoints
    
    dist.type = strsplit(FC, "\\(")[[1]][1]
    dist.number = strsplit(FC, "\\(")[[1]][2]
    dist.number1 = as.numeric(strsplit(dist.number, ",")[[1]][1])
    dist.number2 = as.numeric(strsplit(strsplit(dist.number, 
        ",")[[1]][2], "\\)")[[1]][1])
    
    libsizes = NULL
    mean_fc_relation = NULL
    
    if (is.null(inputCount) == FALSE & is.null(inputLabel) == FALSE) {
        print(paste("Learning parameters from", inputCount))
        model.data <- read.delim(inputCount, header = TRUE, stringsAsFactors = FALSE, 
            row.names = 1)
        
        ind <- apply(model.data, 1, mean) > minAbund
        model.data <- model.data[ind, ]
        
        count.control <- model.data  #[,inputLabel==0]
        ND.control <- apply(count.control, 2, mean)
        count.control <- sweep(count.control, 2, ND.control, 
            "/") * mean(ND.control)
        mu.control <- apply(count.control, 1, mean)
        
        model.input <- cbind(1:nrow(model.data), mu.control)
        model.input <- model.input[order(model.input[, 2], decreasing = TRUE),]
        
        EC_data = dim(model.data)[1]
        ND_data = mean(mu.control)
        
        if (SimulType == "auto") {
            EC = dim(model.data)[1]
            NR1 = length(inputLabel[inputLabel == 0])
            NR2 = length(inputLabel[inputLabel == 1])
            ND = mean(mu.control)
            libsizes = apply(model.data, 2, sum)
        }
        
        NR1= length(inputLabel[inputLabel==0])
        NR2=length(inputLabel[inputLabel==1])
        
        para <- tryCatch(learn_parameter_DESeq(model.data, inputLabel), error=function(e) e)
        
        if (is.null(para$res)){
          para <- learn_parameter_edgeR(model.data, inputLabel)
        }
        
        res <- para$res
        dispersions <- para$dispersions
        
        mean_fc_relation = data.frame(meanA=res$baseMeanA, meanB=res$baseMeanB)
        
        row.names(mean_fc_relation) <- row.names(res)
        
        idx_up <- NULL
        idx_dn <- NULL
        
        ttl_num = 0
        for (i in 1:dim(res)[1]) {
            if(is.na(res$padj[i])) {
              next
            }else if (runif(1) > res$padj[i]){
              ttl_num = ttl_num + 1
            }
        }
        percentFDR <- ttl_num/dim(res)[1]
        
        
        ttl_num = 0
        for (i in 1:dim(res)[1]) {
          if(is.na(res$padj[i])) {
            next
          }else if(runif(1) > res$padj[i] && res$padj[i] < 0.05) {
                ttl_num = ttl_num + 1
          }
        }
        
        FDR5percent <- ttl_num/dim(res)[1]
        
        add_percent <- max(0, percentFDR - FDR5percent)
        add_num <- round(add_percent * (dim(res)[1]))
        
        
        add_FC <- seq(1, dim(res)[1])[abs(res$log2FoldChange) > log2(1.5) & 
            res$padj >= 0.05]
        add_FC <- sample(add_FC, min(length(add_FC), add_num))
        
        for (i in 1:dim(res)[1]) {
          if(is.na(res$padj[i])) {
            next
          }else if ((runif(1) > res$padj[i] && res$padj[i] < 0.05) || 
                i %in% add_FC) {
                if (res$log2FoldChange[i] > 0) {
                  idx_up <- c(idx_up, i)
                }
                if (res$log2FoldChange[i] < 0) {
                  idx_dn <- c(idx_dn, i)
                }
            }
        }
        
        
        DEidx <- c(idx_up, idx_dn)
        NDEidx <- seq(1, dim(res)[1])[-c(idx_up, idx_dn)]
        
        sig.data <- res[c(idx_up, idx_dn), ]
        up.data <- res[idx_up, ]
        dn.data <- res[idx_dn, ]
        
        minFC = min(c(2^(up.data[,"log2FoldChange"]),2^(-dn.data[,"log2FoldChange"])))
        minFC = min(minFC, 1)
        
        fcList=c(2^(up.data[,"log2FoldChange"]),2^(-dn.data[,"log2FoldChange"]))
        
        
        fcList[fcList < minFC] <- minFC
        fcList = sample(fcList, length(fcList))
        
        perDE = dim(sig.data)[1]/dim(model.data)[1]
        upDE = dim(up.data)[1]/dim(model.data)[1]
        downDE = dim(dn.data)[1]/dim(model.data)[1]
        
        
        
        if (SimulType == "auto1") {
            EC = dim(model.data)[1]
            NR1 = ControlRep
            NR2 = CaseRep
            ND = numDataPoints
            DataLib = apply(model.data, 2, sum)
            libsizes <- round(runif(NR1 + NR2, EC * ND * min(DataLib)/mean(DataLib), 
                EC * ND * max(DataLib)/mean(DataLib)))
            
            mean_fc_relation[,1] = mean_fc_relation[,1]*ND/mean(mean_fc_relation[,1])
            mean_fc_relation[,2] = mean_fc_relation[,2]*ND/mean(mean_fc_relation[,2])
        }
        
        if (SimulType == "auto2") {
            EC = EntityCount
            NR1 = ControlRep
            NR2 = CaseRep
            ND = numDataPoints
            DataLib = apply(model.data, 2, sum)
            libsizes <- round(runif(NR1 + NR2, EC * ND * min(DataLib)/mean(DataLib), 
                EC * ND * max(DataLib)/mean(DataLib)))
            
            ratio = EC/dim(mean_fc_relation)[1]
            numDE = round(length(DEidx) * ratio)
            numNDE = EC - numDE
            
            if (ratio > 1) {
                DEidx <- sample(DEidx, numDE, replace = TRUE)
                NDEidx <- sample(NDEidx, numNDE, replace = TRUE)
            } else {
                DEidx <- sample(DEidx, numDE, replace = FALSE)
                NDEidx <- sample(NDEidx, numNDE, replace = FALSE)
            }
            
            fcList <- mean_fc_relation[DEidx, 1]
            mean_fc_relation <- mean_fc_relation[c(DEidx, NDEidx),]
            mean_fc_relation[,1] = mean_fc_relation[,1]*ND/mean(mean_fc_relation[,1])
            mean_fc_relation[,2] = mean_fc_relation[,2]*ND/mean(mean_fc_relation[,2])
            dispersions <- dispersions[c(DEidx, NDEidx),]
            
            DEidx <- 1:length(DEidx)
        }
        
    } else {
        if (dist.type == "Unif") {
            minFC = dist.number1
            maxFC = dist.number2
            fcList = runif(EC * perDE, minFC, maxFC)
        }else if (dist.type == "Norm") {
            fcList = rnorm(EC * perDE, dist.number1, dist.number2)
            minFC = min(fcList)
            maxFC = max(fcList)
        }else if (dist.type == "log2Norm") {
          fcList = 2^(rnorm(EC * perDE, dist.number1, dist.number2))
          minFC = min(fcList)
          maxFC = max(fcList)
        }else if (dist.type == "logNorm") {
          fcList = exp(rnorm(EC * perDE, dist.number1, dist.number2))
          minFC = min(fcList)
          maxFC = max(fcList)
        }
        
        fcList[fcList < 1] <- 1
        
        theta <- 0.5
        if (is.character(SampleVar))
        {
          switch(SampleVar, low = {
            k = 0.05
          }, medium = {
            k = 0.5
          }, high = {
            k = 0.85
          })
        } else {
            k = SampleVar * 2
        }
        
        dispersions <- rgamma(EC, shape = k, scale = theta)
        dispersions <- cbind(dispersions, dispersions)
        
        if (AbundProfile == "HBR") {
            data("HBR")
            model.input = HBR
            print(paste("Using", AbundProfile, "Profile"))
        } else if (AbundProfile == "BP")
        {
          data("BP")
          model.input = BP
          print(paste("Using", AbundProfile, "Profile"))
        } else if (AbundProfile == "Wu")
        {
          data("Wu")
          model.input = Wu
          print(paste("Using", AbundProfile, "Profile"))
        }else{
          
          if (AbundProfile == "SingleCell"){
            data(SingleCell)
            model.data <- SingleCell
          }else{
            model.data <- read.delim(AbundProfile, header = TRUE, stringsAsFactors = FALSE, 
                                     row.names = 1)
          }

          if (dim(model.data)[2]>1)
          {
            ind <- apply(model.data, 1, mean) > minAbund
            model.data <- model.data[ind, ]
            
            count.control <- model.data  
            ND.control <- apply(count.control, 2, mean)
            count.control <- sweep(count.control, 2, ND.control, 
                                   "/") * mean(ND.control)
            mu.control <- apply(count.control, 1, mean)
            
            print(paste("Using Abundance Profile obtained from", AbundProfile))
            
          }else{
            ind <- model.data[,1] > minAbund
            model.data <- as.data.frame(model.data[ind, ])
            mu.control <- model.data[,1]
            print(paste("Using", AbundProfile, "Profile"))
          }

          
          model.input <- cbind(1:nrow(model.data), mu.control)
          model.input <- model.input[order(model.input[, 2], decreasing = TRUE),]
         
        }

    }
    
    dataLabel = c(rep(0, NR1), rep(1, NR2))
    
    # Step 2: Randomly select genes
    geneCounts <- EC
    lowerLim <- round(geneCounts * 0.1)
    upperLimL <- tail(which(model.input[, 2] > ((1 * minAbund) - 
        1)), n = 1) - round(geneCounts * 0.1) + 1
    upperLimR <- tail(which(model.input[, 2] > ((1 * minAbund) - 
        1)), n = 1)
    
    # Selecting top lowerLim expressed genes
    model.selected <- model.input[1:lowerLim, ]
    
    # Selecting bottom lowerLim expressed genes
    model.selected <- rbind(model.selected, model.input[upperLimL:upperLimR, 
        ])
    
    # Randomly selecting remaining genes
    if (length((lowerLim + 1):(upperLimL - 1)) >= (geneCounts - 
        2 * lowerLim)) {
        randomList <- sample((lowerLim + 1):(upperLimL - 1), 
            (geneCounts - 2 * lowerLim), replace = FALSE)
    } else {
        randomList <- sample((lowerLim + 1):(upperLimL - 1), 
            (geneCounts - 2 * lowerLim), replace = TRUE)
    }
    model.selected <- rbind(model.selected, model.input[randomList, 
        ])
    rownames(model.selected) <- paste("g", 1:dim(model.selected)[1], 
        sep = "")
    # rownames(model.selected) <-model.selected[,1];
    
    model.temp <- as.matrix(model.selected[, -1])
    model.rawFreq <- sweep(model.temp, 2, apply(model.temp, 2, 
        sum), "/")
    
    model.rawFreq <- sample(model.rawFreq, dim(model.rawFreq)[1])
    model.matrix <- matrix(0, nrow = EC, ncol = (NR1 + NR2))
    rownames(model.matrix) <- paste("g", 1:dim(model.selected)[1], 
        sep = "")
    
    
    # DE simulation
    EDgenelist <- data.frame(matrix(0, nrow = length(fcList), 
        ncol = 4))
    colnames(EDgenelist) <- c("RowID", "geneName", "FC", "log2FC")
    
    randomList <- sample(1:dim(model.matrix)[1], length(fcList), 
        replace = FALSE)
    coinList <- c(rep(1, ceiling(EC * upDE) + 1), rep(0, (EC * 
        downDE) + 1))
    coinList <- sample(coinList, length(fcList))
    # coinList <- as.matrix(transform(sample(coinList)));
    index <- 1
    
    
    if (is.null(libsizes)) 
        libsizes <- round(runif(NR1 + NR2, EC * ND * (1 - varLibsizes), 
            EC * ND * (1 + varLibsizes)))
    
    if (is.null(mean_fc_relation))
    {
      AllFC <- rep(1, length(model.rawFreq)) 
      fcList[coinList==0] <- 1/fcList[coinList==0]
      AllFC[randomList] <- fcList
      
      meanA <- NULL
      meanB <- NULL
      for(j in 1:length(model.rawFreq)){
        mean12 <- model.rawFreq[j]*mean(libsizes)
        fc <- AllFC[j]
        if (fc>=1)
        {
          coin <- 0
        }else if(fc<1)
        {
          coin <- 1
          fc <- 1/fc
        }
        mu1 = (coin * sqrt(fc) + (1-coin) / sqrt(fc)) * mean12
        mu2 = ((1-coin) * sqrt(fc) + coin / sqrt(fc)) * mean12
        meanA <- c(meanA, mu1)
        meanB <- c(meanB, mu2)
      }
      
      mean_fc_relation = data.frame(meanA=meanA, meanB=meanB)
      DEidx <- randomList
    }
    
    if (model == "NegBinomial" || model == "Full") {
        print(paste("Simulation model: ", model))
        for (j in 1:EC) {
            if (j %in% DEidx) {
              mean1 <- mean_fc_relation[j,1]
              mean2 <- mean_fc_relation[j,2]
              fc <- mean2/mean1
                
              for(i in 1:NR1)
                model.matrix[j,i] <- rnbinom(1, size = 1/dispersions[j,1], 
                                             mu = mean1)      
              for(i in (NR1+1):(NR1+NR2))
                model.matrix[j,i] <- rnbinom(1, size = 1/dispersions[j,2], 
                                             mu = mean2)
              
              EDgenelist[index,1]=j
              EDgenelist[index,2]=rownames(model.matrix)[j]
              EDgenelist[index,3]=fc
              EDgenelist[index,4]=log2(fc)
              
              index <- index + 1
            } else {
                
              mean1 <- mean_fc_relation[j,1]
              mean2 <- mean_fc_relation[j,2]
              
              for(i in 1:NR1)
                model.matrix[j,i] <- rnbinom(1, size = 1/dispersions[j,1], 
                                             mu = (mean1+mean2)/2);      
              for(i in (NR1+1):(NR1+NR2))
                model.matrix[j,i] <- rnbinom(1, size = 1/dispersions[j,2], 
                                             mu = (mean1+mean2)/2)
            }
            
        }
    } else if (model == "Multinom") {
        print("Simulation model: Multinomial")

        for (j in 1:EC) {
            if (j %in% DEidx) 
            {
                mean1 <- mean_fc_relation[j,1]
                mean2 <- mean_fc_relation[j,2]
                fc <- mean2/mean1
                
                for (i in 1:NR1) 
                {
                  model.matrix[j, i] <- mean1/mean(libsizes)*libsizes[i]
                }
                    
                for (i in (NR1 + 1):(NR1 + NR2)) {
                    model.matrix[j, i] <-  mean2/mean(libsizes)*libsizes[i]
                }
                EDgenelist[index,1]=j
                EDgenelist[index,2]=rownames(model.matrix)[j]
                EDgenelist[index,3]=fc
                EDgenelist[index,4]=log2(fc)
                
                index <- index + 1
            } else {
              mean1 <- mean_fc_relation[j,1]
              mean2 <- mean_fc_relation[j,2]
              
              mean12 <- (mean1+mean2)/2
              for (i in 1:NR1) 
              {
                model.matrix[j, i] <- mean12/mean(libsizes)*libsizes[i]
              }
              
              for (i in (NR1 + 1):(NR1 + NR2)) {
                model.matrix[j, i] <- mean12/mean(libsizes)*libsizes[i]
              }
            }
            
        }
        
    } else if (model == "ModelFree"||SimulModel=="ModelFreeMn") {
        print("Simulation model: Model Free")
        
        if (!is.null(modelFile))
        {
          print(paste("Using", modelFile, "as Model data for sampling"))
          if (modelFile %in% c("SingleCell")){
            data(SingleCell)
            model.data <- SingleCell
          }else if (file.exists(modelFile)){
            model.data <- read.delim(modelFile, header = TRUE, stringsAsFactors = FALSE, 
                                     row.names = 1)
            if (dim(model.data)[2]<=1)
            {
              stop("ERROR: No replicates!")
            }
          }else{
            stop("ERROR: modelFile not correct! Set modelFile to NULL or correct name or correct path!")
          }
        }
        
        if (!is.null(inputCount) && !is.null(inputLabel) && is.null(modelFile))
        {
          cond1 <- 0
          cond2 <- 1
        }else{
          inputLabel <- rep(0, times = dim(model.data)[2])
          cond1 <- 0
          cond2 <- 0
        }

        ind <- apply(model.data, 1, mean) > 0
        model.data <- model.data[ind, ]
        
        cds <- newCountDataSet(model.data, inputLabel)
        cds <- estimateSizeFactors(cds)
        data_norm <- counts(cds, normalized = TRUE)
          
        mean_data1 <- apply(data_norm[, inputLabel == cond1], 1, 
                            mean)
        mean_data2 <- apply(data_norm[, inputLabel == cond2], 1, 
                            mean)
        
        rindex <- 1:dim(data_norm)[1]
        
        ModelNR1 <- sum(inputLabel == cond1)
        ModelNR2 <- sum(inputLabel == cond2)
        
        for (j in 1:EC) {
            
            if (j %in% DEidx) {
              
              mean1 <- mean_fc_relation[j,1]
              mean2 <- mean_fc_relation[j,2]
              fc <- mean2/mean1
              
              tmp <- abs(mean_data1 - mean1)
              if (NR1 < ModelNR1)
              {
                idx <- order(tmp)[1]
                idx1 <- rindex[tmp <= tmp[idx]]
                if (length(idx1)>1){
                  idx<- sample(idx1,1)
                }
                
                d1 <- as.vector(data_norm[idx, inputLabel == cond1])
                if (mean(d1)!=0) d1 <- mean1*d1/mean(d1)
              }else{
                thrsh = min(10, mean1 * 0.1)
                if (sum(tmp < thrsh) > 10) {
                  idx <- rindex[tmp < thrsh]
                } else {
                  idx <- order(tmp)[1:10]
                }
              
                d1 <- data_norm[idx, inputLabel == cond1]
                if (any(apply(d1, 1, mean) == 0)) {
                  d1 <- as.vector(d1)
                } else {
                  d1 <- as.vector(mean1 * d1/apply(d1, 1, mean))

                }
              
              }
              model.matrix[j, 1:NR1] <- sample(d1, NR1, replace = TRUE)
                  
              tmp <- abs(mean_data2 - mean2)
              if (NR2 < ModelNR2)
              {
                idx <- order(tmp)[1]
                idx1 <- rindex[tmp <= tmp[idx]]
                if (length(idx1)>1){
                  idx<- sample(idx1,1)
                }
                d1 <- as.vector(data_norm[idx, inputLabel == cond2])
                if (mean(d1)!=0) d1 <- mean2*d1/mean(d1)
              }else{
                thrsh = min(10, mean2 * 0.1)
                if (sum(tmp < thrsh) > 10) {
                  idx <- rindex[tmp < thrsh]
                } else {
                  idx <- order(tmp)[1:10]
                }
                
                d1 <- data_norm[idx, inputLabel == cond2]
                if (any(apply(d1, 1, mean) == 0)) {
                  d1 <- as.vector(d1)
                } else {
                  d1 <- as.vector(mean2 * d1/apply(d1, 1, mean))
                  
                }
                
              }
              model.matrix[j, (NR1 + 1):(NR1 + NR2)] <- sample(d1,  NR2, replace = TRUE)
              
                
              EDgenelist[index,1]=j
              EDgenelist[index,2]=rownames(model.matrix)[j]
              EDgenelist[index,3]=fc
              EDgenelist[index,4]=log2(fc)
              
              index <- index + 1
            
            } else {
                
              mean1 <- mean_fc_relation[j,1]
              mean2 <- mean_fc_relation[j,2]
              mean12 <- (mean1+mean2)/2
              
              tmp <- abs(mean_data1 - mean12)
              if (NR1 < ModelNR1)
              {
                idx <- order(tmp)[1]
                idx1 <- rindex[tmp <= tmp[idx]]
                if (length(idx1)>1){
                  idx<- sample(idx1,1)
                }
                d1 <- as.vector(data_norm[idx, inputLabel == cond1])
                if (mean(d1)!=0) d1 <- mean12*d1/mean(d1)
              }else{
                thrsh = min(10, mean12 * 0.1)
                if (sum(tmp < thrsh) > 10) {
                  idx <- rindex[tmp < thrsh]
                } else {
                  idx <- order(tmp)[1:10]
                }
                
                d1 <- data_norm[idx, inputLabel == cond1]
                if (any(apply(d1, 1, mean) == 0)) {
                  d1 <- as.vector(d1)
                } else {
                  d1 <- as.vector(mean12 * d1/apply(d1, 1, mean))
                  
                }
                
              }
              model.matrix[j, 1:NR1] <- sample(d1, NR1, replace = TRUE)
              
              tmp <- abs(mean_data2 - mean12)
              if (NR2 < ModelNR2)
              {
                idx <- order(tmp)[1]
                idx1 <- rindex[tmp <= tmp[idx]]
                if (length(idx1)>1){
                  idx<- sample(idx1,1)
                }
                d1 <- as.vector(data_norm[idx, inputLabel == cond2])
                if (mean(d1)!=0) d1 <- mean12*d1/mean(d1)
              }else{
                thrsh = min(10, mean12 * 0.1)
                if (sum(tmp < thrsh) > 10) {
                  idx <- rindex[tmp < thrsh]
                } else {
                  idx <- order(tmp)[1:10]
                }
                
                d1 <- data_norm[idx, inputLabel == cond2]
                if (any(apply(d1, 1, mean) == 0)) {
                  d1 <- as.vector(d1)
                } else {
                  d1 <- as.vector(mean12 * d1/apply(d1, 1, mean))
                  
                }
                
              }
              model.matrix[j, (NR1 + 1):(NR1 + NR2)] <- sample(d1,  NR2, replace = TRUE)
              
                
                
            }
            
            
        }
        
    } 
        
    # Generate frequency matrix
    model.freq <- sweep(model.matrix, 2, apply(model.matrix, 
        2, sum), "/")
    
    if (model == "NegBinomial"|| model == "ModelFree") 
        model.combined <- round(sweep(model.freq, 2, libsizes, 
            "*"))
    
    if (model == "Multinom" || model == "Full"|| model == "ModelFreeMn") {
        # Generate multinomial distribution using frequen matrix
        model.combined <- matrix(0, nrow = EC, ncol = (NR1 + 
            NR2))
        rownames(model.combined) <- rownames(model.matrix)
        for (i in 1:(NR1 + NR2)) {
            model.combined[, i] <- rmultinom(n = 1, size = libsizes[i], 
                prob = model.freq[, i])
        }
    }
    
    
    if (outlier==TRUE && SimulModel != "ModelFree" && SimulModel != "ModelFreeMn"){
      numG=dim(model.combined)[1]
      numR=dim(model.combined)[2]
      ind <- sample(1:numG, round(numG*perOutlier))
      for (i in ind){
        j <- sample(numR,1)
        if (runif(1)>0.5) {
          model.combined[i,j] <- round(model.combined[i,j]*factorOutlier)
        }else{
          model.combined[i,j] <- round(model.combined[i,j]/factorOutlier)
        }
        
      }
    }
    return(list(count = model.combined, DiffAbundList = EDgenelist, 
        dataLabel = dataLabel))
}


learn_parameter_DESeq <- function(data, inputLabel){
  NR1= length(inputLabel[inputLabel==0])
  NR2=length(inputLabel[inputLabel==1])
  
  runID=""
  
  conds <- c(rep("N",NR1), rep("T", NR2));
  cds <- newCountDataSet(normalizeData(data, conds, runID, 10)$normCounts, conditions=inputLabel)
  sizeFactors(cds) <- c(rep(1,length(conds)))
  #cds <- newCountDataSet(model.data, conditions=inputLabel)
  #cds <- estimateSizeFactors(cds)
  
  
  if(length(inputLabel) == 2){
    cds <- estimateDispersions(cds, method= "blind", fitType='local', sharingMode='fit-only')
  }else{
    cds <- estimateDispersions(cds, method= "per-condition", fitType='local')
  }
  
  ## differential expression
  res <- nbinomTest(cds, "0", "1")
  rownames(res) <- res$id
  
#  mean_fc_relation = data.frame(meanA=res$baseMeanA, meanB=res$baseMeanB)
#  row.names(mean_fc_relation) <- row.names(res)
  
  dispersions <- fData(cds)
  
  return(list(res=res, dispersions=dispersions))
}


learn_parameter_edgeR <- function(data, inputLabel){

  NR1= length(inputLabel[inputLabel==0])
  NR2=length(inputLabel[inputLabel==1])
  
  runID=""
  
  conds <- c(rep("N",NR1), rep("T", NR2));
  model.data.norm <- normalizeData(data, conds, runID, 10)$normCounts
  
  d <- DGEList(counts=model.data.norm, group=inputLabel)
  d$samples$norm.factors <- rep(1, length(inputLabel));
  
  if(2 == length(inputLabel)){
    d$common.dispersion = 0.4
    
  }else{
    ## estimate common dispersion
    d <- estimateCommonDisp(d)
    
    ## estimate gene specific dispersion
    d <- estimateTagwiseDisp(d)
  }
  
  de.com = exactTest(d);
  res <- topTags(de.com, n=nrow(de.com))$table;
  
  res$padj <- res$FDR
  res$log2FoldChange <- res$logFC
  
  
  res$baseMeanA = apply(as.matrix(model.data.norm[,inputLabel==0]),1,mean)[row.names(res)]
  res$baseMeanB = apply(as.matrix(model.data.norm[,inputLabel==1]),1,mean)[row.names(res)]
  
  dispersions <- cbind(d$tagwise.dispersion,d$tagwise.dispersion)
  rownames(dispersions) <- rownames(d$counts)
  dispersions <- dispersions[rownames(mean_fc_relation),]
  
  return(list(res=res, dispersions=dispersions))
}