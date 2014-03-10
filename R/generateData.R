generateData <- function(SimulModel = "Full", SampleVar = "medium", 
    ControlRep = 5, CaseRep = ControlRep, EntityCount = 1000, 
    FC = "Norm(2,1)", perDiffAbund = 0.1, upPDA = perDiffAbund/2, 
    downPDA = perDiffAbund/2, numDataPoints = 100, modelFile = "HBRmodel", 
    minAbund = 10, varLibsizes = 0.1, inputCount = NULL, inputLabel = NULL, 
    SimulType = "auto") {
    
    NDEfc = FALSE
    
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
        print(paste("Using", inputCount))
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
        model.input <- model.input[order(model.input[, 2], decreasing = TRUE), 
            ]
        
        EC_data = dim(model.data)[1]
        ND_data = mean(mu.control)
        
        if (SimulType == "auto") {
            EC = dim(model.data)[1]
            NR1 = length(inputLabel[inputLabel == 0])
            NR2 = length(inputLabel[inputLabel == 1])
            ND = mean(mu.control)
            libsizes = apply(model.data, 2, sum)
        }
        
        # library('edgeR')
        
        d <- DGEList(counts = model.data, group = inputLabel)
        d = calcNormFactors(d, method = "TMM")
        
        if (2 == length(inputLabel)) {
            d$common.dispersion = 0.4
            
        } else {
            ## estimate common dispersion
            d <- estimateCommonDisp(d)
            
            ## estimate gene specific dispersion
            d <- estimateTagwiseDisp(d)
        }
        
        de.com = exactTest(d)
        res <- topTags(de.com, n = nrow(de.com))$table
        
        mean_fc_relation = data.frame(mean = 2^res$logCPM, FC = 2^res$logFC)
        
        ### 
        mean_fc_relation[, 1] = mean_fc_relation[, 1] * EC_data * 
            ND_data/10^6
        ### 
        row.names(mean_fc_relation) <- row.names(res)
        
        dispersions <- d$tagwise.dispersion
        names(dispersions) <- rownames(d$counts)
        dispersions <- dispersions[rownames(mean_fc_relation)]
        
        idx_up <- NULL
        idx_dn <- NULL
        
        ttl_num = 0
        for (i in 1:dim(res)[1]) {
            if (runif(1) > res$FDR[i]) 
                ttl_num = ttl_num + 1
        }
        percentFDR <- ttl_num/dim(res)[1]
        
        
        ttl_num = 0
        for (i in 1:dim(res)[1]) {
            if (runif(1) > res$FDR[i] && res$FDR[i] < 0.05) {
                ttl_num = ttl_num + 1
            }
        }
        
        FDR5percent <- ttl_num/dim(res)[1]
        
        add_percent <- max(0, percentFDR - FDR5percent)
        add_num <- round(add_percent * (dim(res)[1]))
        
        
        add_FC <- seq(1, dim(res)[1])[abs(res$logFC) > log2(1.5) & 
            res$FDR >= 0.05]
        add_FC <- sample(add_FC, min(length(add_FC), add_num))
        
        for (i in 1:dim(res)[1]) {
            if ((runif(1) > res$FDR[i] && res$FDR[i] < 0.05) || 
                i %in% add_FC) {
                if (res$logFC[i] > 0) {
                  idx_up <- c(idx_up, i)
                }
                if (res$logFC[i] < 0) {
                  idx_dn <- c(idx_dn, i)
                }
            }
        }
        
        
        DEidx <- c(idx_up, idx_dn)
        NDEidx <- seq(1, dim(res)[1])[-c(idx_up, idx_dn)]
        nonDEfc <- res$logFC[-c(idx_up, idx_dn)]
        
        # print(length(idx_up)) print(length(idx_dn))
        # print(length(nonDEfc))
        
        sig.data <- res[c(idx_up, idx_dn), ]
        up.data <- res[idx_up, ]
        dn.data <- res[idx_dn, ]
        
        # print(2^min(abs(sig.data$logFC)))
        
        minFC = min(c(2^(up.data[, 1]), 2^(-dn.data[, 1])))
        minFC = min(minFC, 1)
        
        fcList = c(2^(up.data[, 1]), 2^(-dn.data[, 1]))
        
        
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
            # print(libsizes)
            mean_fc_relation[, 1] = mean_fc_relation[, 1] * ND/mean(mean_fc_relation[, 
                1])
        }
        
        if (SimulType == "auto2") {
            EC = EntityCount
            NR1 = ControlRep
            NR2 = CaseRep
            ND = numDataPoints
            DataLib = apply(model.data, 2, sum)
            libsizes <- round(runif(NR1 + NR2, EC * ND * min(DataLib)/mean(DataLib), 
                EC * ND * max(DataLib)/mean(DataLib)))
            # print(libsizes)
            
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
            
            fcList <- mean_fc_relation[DEidx, 2]
            mean_fc_relation <- mean_fc_relation[c(DEidx, NDEidx), 
                ]
            mean_fc_relation[, 1] = mean_fc_relation[, 1] * ND/mean(mean_fc_relation[, 
                1])
            dispersions <- dispersions[c(DEidx, NDEidx)]
            
            DEidx <- 1:length(DEidx)
        }
        
    } else {
        if (dist.type == "Unif") {
            minFC = dist.number1
            maxFC = dist.number2
            fcList = runif(EC * perDE, minFC, maxFC)
        }
        if (dist.type == "Norm") {
            fcList = 2^(rnorm(EC * perDE, dist.number1, dist.number2))
            minFC = min(fcList)
            maxFC = max(fcList)
        }
        
        fcList[fcList < 1] <- 1
        nonDEfc <- c(0, 0)
        
        theta <- 0.5
        if (is.character(SampleVar)) 
            switch(SampleVar, low = {
                k = 0.05
            }, medium = {
                k = 0.5
            }, high = {
                k = 0.85
            }) else {
            k = SampleVar * 2
        }
        
        dispersions <- rgamma(EC, shape = k, scale = theta)
        if (modelFile == "HBRmodel") {
            #HBRmodel <- read.table(system.file("extdata","HBRmodel.txt.gz", 
             #           package="EDDA"),head=T)
            #rm(HBRmodel)
            data("HBRmodel")
            #load(system.file("extdata","HBRmodel.txt.gz", package="EDDA")) 
            model.input = HBRmodel
            print(paste("Using", modelFile))
        } else model.input <- read.delim(modelFile, header = FALSE, 
            stringsAsFactors = TRUE)
    }
    
    dataLabel = c(rep(0, NR1), rep(1, NR2))
    y <- matrix(0, nrow = EC, ncol = (NR1 + NR2))
    
    
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
    DEgeneID <- NULL
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
    
    if (is.null(mean_fc_relation)) {
        AllFC <- rep(1, length(model.rawFreq))
        fcList[coinList == 0] <- 1/fcList[coinList == 0]
        AllFC[randomList] <- fcList
        mean_fc_relation = data.frame(mean = model.rawFreq * 
            mean(libsizes), FC = AllFC)
        DEidx <- randomList
    }
    
    if (model == "NegBinomial" || model == "Full") {
        print(paste("Simulation model: ", model))
        for (j in 1:EC) {
            if (j %in% DEidx) {
                mean12 <- mean_fc_relation[j, 1]
                fc <- mean_fc_relation[j, 2]
                if (fc >= 1) {
                  coin <- 0
                } else if (fc < 1) {
                  coin <- 1
                  fc <- 1/fc
                }
                
                # coin <- coinList[index] # 0: up-regulation; 1:
                # down-regulation fc <- fcList[index]
                
                for (i in 1:NR1) model.matrix[j, i] <- rnbinom(1, 
                  size = 1/dispersions[j], mu = (coin * sqrt(fc) + 
                    (1 - coin)/sqrt(fc)) * mean12)
                for (i in (NR1 + 1):(NR1 + NR2)) model.matrix[j, 
                  i] <- rnbinom(1, size = 1/dispersions[j], mu = ((1 - 
                  coin) * sqrt(fc) + coin/sqrt(fc)) * mean12)
                if (coin == 1) {
                  EDgenelist[index, 1] = j
                  EDgenelist[index, 2] = rownames(model.matrix)[j]
                  EDgenelist[index, 3] = fc
                  EDgenelist[index, 4] = log2(fc)
                } else {
                  EDgenelist[index, 1] = j
                  EDgenelist[index, 2] = rownames(model.matrix)[j]
                  EDgenelist[index, 3] = 1/fc
                  EDgenelist[index, 4] = log2(1/fc)
                }
                DEgeneID <- c(DEgeneID, rownames(model.matrix)[j])
                index <- index + 1
            } else {
                
                mean12 <- mean_fc_relation[j, 1]
                if (NDEfc == TRUE) {
                  fc <- mean_fc_relation[j, 2]
                } else {
                  fc = 1
                }
                # coin = sample(c(0,1),1)
                if (fc >= 1) {
                  coin <- 0
                } else if (fc < 1) {
                  coin <- 1
                  fc <- 1/fc
                }
                
                
                for (i in 1:NR1) model.matrix[j, i] <- rnbinom(1, 
                  size = 1/dispersions[j], mu = (coin * sqrt(fc) + 
                    (1 - coin)/sqrt(fc)) * mean12)
                for (i in (NR1 + 1):(NR1 + NR2)) model.matrix[j, 
                  i] <- rnbinom(1, size = 1/dispersions[j], mu = ((1 - 
                  coin) * sqrt(fc) + coin/sqrt(fc)) * mean12)
            }
            
        }
    } else if (model == "Multinom") {
        print("Simulation model: Multinomial")
        
        for (j in 1:EC) {
            if (j %in% randomList) 
                {
                  coin <- coinList[index]  # 0: up-regulation; 1: down-regulation
                  fc <- fcList[index]
                  for (i in 1:NR1) {
                    model.matrix[j, i] <- (coin * sqrt(fc) + 
                      (1 - coin)/sqrt(fc)) * model.rawFreq[j] * 
                      libsizes[i]
                  }
                  for (i in (NR1 + 1):(NR1 + NR2)) {
                    model.matrix[j, i] <- ((1 - coin) * sqrt(fc) + 
                      coin/sqrt(fc)) * model.rawFreq[j] * libsizes[i]
                  }
                  if (coin == 1) {
                    EDgenelist[index, 1] = j
                    EDgenelist[index, 2] = rownames(model.matrix)[j]
                    EDgenelist[index, 3] = fc
                    EDgenelist[index, 4] = log2(fc)
                  } else {
                    EDgenelist[index, 1] = j
                    EDgenelist[index, 2] = rownames(model.matrix)[j]
                    EDgenelist[index, 3] = 1/fc
                    EDgenelist[index, 4] = log2(1/fc)
                  }
                  DEgeneID <- c(DEgeneID, rownames(model.matrix)[j])
                  index <- index + 1
                }  # end for loop
 else {
                for (i in 1:(NR1 + NR2)) model.matrix[j, i] <- model.rawFreq[j] * 
                  libsizes[i]
            }
            
        }
        
        
    } else if (model == "ModelFree") {
        print("Simulation model: Model Free")
        # library(DESeq)
        cds <- newCountDataSet(model.data, inputLabel)
        cds <- estimateSizeFactors(cds)
        data_norm <- round(counts(cds, normalized = TRUE))
        
        mean_data0 <- apply(data_norm[, inputLabel == 0], 1, 
            mean)
        mean_data1 <- apply(data_norm[, inputLabel == 1], 1, 
            mean)
        
        data0_libsize <- mean(apply(data_norm[, inputLabel == 
            0], 2, sum))
        data1_libsize <- mean(apply(data_norm[, inputLabel == 
            1], 2, sum))
        
        rindex <- 1:dim(data_norm)[1]
        
        for (j in 1:EC) {
            
            if (j %in% DEidx) {
                mean12 <- mean_fc_relation[j, 1]
                fc <- mean_fc_relation[j, 2]
                if (fc >= 1) {
                  coin <- 0
                } else if (fc < 1) {
                  coin <- 1
                  fc <- 1/fc
                }
                # coin <- coinList[index] # 0: up-regulation; 1:
                # down-regulation fc <- fcList[index]
                
                {
                  mu = (coin * sqrt(fc) + (1 - coin)/sqrt(fc)) * 
                    mean12
                  # mu = mean_data0[row.names(mean_fc_relation)[j]] mu =
                  # mean_data0[paste('g',j,sep='')]
                  
                  tmp <- abs(mean_data0 - mu)
                  thrsh = min(10, mu * 0.1)
                  if (sum(tmp < thrsh) > 10) {
                    idx <- rindex[tmp < thrsh]
                  } else {
                    idx <- order(tmp)[1:10]
                  }
                  
                  d1 <- data_norm[idx, inputLabel == 0]
                  if (any(apply(d1, 1, mean) == 0)) {
                    d1 <- as.vector(d1)
                  } else {
                    d1 <- round(as.vector(mu * d1/apply(d1, 1, 
                      mean)))
                    # d1 <- as.vector(d1)
                  }
                  model.matrix[j, 1:NR1] <- sample(d1, NR1, replace = TRUE)
                  
                }
                
                
                {
                  mu = ((1 - coin) * sqrt(fc) + coin/sqrt(fc)) * 
                    mean12
                  # mu = mean_data1[row.names(mean_fc_relation)[j]] mu =
                  # mean_data1[paste('g',j,sep='')]
                  
                  tmp <- abs(mean_data1 - mu)
                  thrsh = min(10, mu * 0.1)
                  if (sum(tmp < thrsh) > 10) {
                    idx <- rindex[tmp < thrsh]
                  } else {
                    idx <- order(tmp)[1:10]
                  }
                  
                  d1 <- data_norm[idx, inputLabel == 1]
                  if (any(apply(d1, 1, mean) == 0)) {
                    d1 <- as.vector(d1)
                  } else {
                    d1 <- round(as.vector(mu * d1/apply(d1, 1, 
                      mean)))
                    # d1 <- as.vector(d1)
                  }
                  model.matrix[j, (NR1 + 1):(NR1 + NR2)] <- sample(d1, 
                    NR2, replace = TRUE)
                  
                }
                
                
                if (coin == 1) {
                  EDgenelist[index, 1] = j
                  EDgenelist[index, 2] = rownames(model.matrix)[j]
                  EDgenelist[index, 3] = fc
                  EDgenelist[index, 4] = log2(fc)
                } else {
                  EDgenelist[index, 1] = j
                  EDgenelist[index, 2] = rownames(model.matrix)[j]
                  EDgenelist[index, 3] = 1/fc
                  EDgenelist[index, 4] = log2(1/fc)
                }
                DEgeneID <- c(DEgeneID, rownames(model.matrix)[j])
                index <- index + 1
                
            } else {
                
                mean12 <- mean_fc_relation[j, 1]
                if (NDEfc == TRUE) {
                  fc <- mean_fc_relation[j, 2]
                } else {
                  fc = 1
                }
                
                if (fc >= 1) {
                  coin <- 0
                } else if (fc < 1) {
                  coin <- 1
                  fc <- 1/fc
                }
                
                {
                  mu = (coin * sqrt(fc) + (1 - coin)/sqrt(fc)) * 
                    mean12
                  # mu = mean_data0[row.names(mean_fc_relation)[j]] mu =
                  # mean_data0[paste('g',j,sep='')]
                  
                  
                  tmp <- abs(mean_data0 - mu)
                  thrsh = min(10, mu * 0.1)
                  if (sum(tmp < thrsh) > 10) {
                    idx <- rindex[tmp < thrsh]
                  } else {
                    idx <- order(tmp)[1:10]
                  }
                  
                  d1 <- data_norm[idx, inputLabel == 0]
                  if (any(apply(d1, 1, mean) == 0)) {
                    d1 <- as.vector(d1)
                  } else {
                    d1 <- round(as.vector(mu * d1/apply(d1, 1, 
                      mean)))
                    # d1 <- as.vector(d1)
                  }
                  model.matrix[j, 1:NR1] <- sample(d1, NR1, replace = TRUE)
                  
                }
                
                
                {
                  mu = ((1 - coin) * sqrt(fc) + coin/sqrt(fc)) * 
                    mean12
                  # mu = mean_data1[row.names(mean_fc_relation)[j]]
                  
                  # mu = mean_data1[paste('g',j,sep='')]
                  
                  tmp <- abs(mean_data1 - mu)
                  thrsh = min(10, mu * 0.1)
                  if (sum(tmp < thrsh) > 10) {
                    idx <- rindex[tmp < thrsh]
                  } else {
                    idx <- order(tmp)[1:10]
                  }
                  
                  d1 <- data_norm[idx, inputLabel == 1]
                  if (any(apply(d1, 1, mean) == 0)) {
                    d1 <- as.vector(d1)
                  } else {
                    d1 <- round(as.vector(mu * d1/apply(d1, 1, 
                      mean)))
                    # d1 <- as.vector(d1)
                  }
                  model.matrix[j, (NR1 + 1):(NR1 + NR2)] <- sample(d1, 
                    NR2, replace = TRUE)
                  
                }
                
            }
            
            
        }
        
    } else if (model == "SingleCell") {
        print("Simulation model: Single Cell")
        # library(DESeq)
        #SingleCell<- read.table(system.file("extdata","SingleCell.txt.gz", 
        #                package="EDDA"),head=T)
        #rm(SingleCell)
        data(SingleCell)
        data <- SingleCell
        
        cds <- newCountDataSet(data, rep(0, times = dim(data)[2]))
        cds <- estimateSizeFactors(cds)
        data_norm <- round(counts(cds, normalized = TRUE))
        
        ind <- apply(data_norm, 1, mean) > 0
        data_norm <- data_norm[ind, ]
        
        mean_ref <- apply(data_norm, 1, mean)
        
        rindex <- 1:dim(data_norm)[1]
        
        for (j in 1:EC) {
            
            if (j %in% DEidx) {
                mean12 <- mean_fc_relation[j, 1]
                fc <- mean_fc_relation[j, 2]
                if (fc >= 1) {
                  coin <- 0
                } else if (fc < 1) {
                  coin <- 1
                  fc <- 1/fc
                }
                # coin <- coinList[index] # 0: up-regulation; 1:
                # down-regulation fc <- fcList[index]
                
                {
                  mu = (coin * sqrt(fc) + (1 - coin)/sqrt(fc)) * 
                    mean12
                  tmp <- abs(mean_ref - mu)
                  idx <- order(tmp)[1]
                  
                  d1 <- as.numeric(data_norm[idx, ])
                  d1 <- round(d1 * mu/mean(d1))
                  
                  stepp = as.numeric(quantile(d1, probs = seq(0, 
                    1, length.out = (NR1 + 1))))
                  
                  for (nr in 1:NR1) {
                    if (nr == 1) {
                      ind <- d1 >= stepp[nr] & d1 <= stepp[nr + 
                        1]
                    } else {
                      ind <- d1 > stepp[nr] & d1 <= stepp[nr + 
                        1]
                    }
                    
                    if (sum(ind) > 0) 
                      model.matrix[j, nr] <- sample(d1[ind], 
                        1, replace = TRUE) else {
                      model.matrix[j, nr] <- stepp[nr]
                    }
                  }
                  
                  d1 <- as.numeric(model.matrix[j, 1:NR1])
                  model.matrix[j, 1:NR1] <- sample(round(d1), 
                    NR1)
                  
                  
                }
                
                
                {
                  mu = ((1 - coin) * sqrt(fc) + coin/sqrt(fc)) * 
                    mean12
                  
                  tmp <- abs(mean_ref - mu)
                  idx <- order(tmp)[1]
                  
                  d1 <- as.numeric(data_norm[idx, ])
                  d1 <- round(d1 * mu/mean(d1))
                  
                  # model.matrix[j,(NR1+1):(NR1+NR2)] <- sample(d1, NR2,
                  # replace=T)
                  
                  stepp = as.numeric(quantile(d1, probs = seq(0, 
                    1, length.out = (NR2 + 1))))
                  
                  for (nr in 1:NR2) {
                    
                    if (nr == 1) {
                      ind <- d1 >= stepp[nr] & d1 <= stepp[nr + 
                        1]
                    } else {
                      ind <- d1 > stepp[nr] & d1 <= stepp[nr + 
                        1]
                    }
                    
                    if (sum(ind) > 0) 
                      model.matrix[j, nr + NR1] <- sample(d1[ind], 
                        1, replace = TRUE) else {
                      model.matrix[j, nr + NR1] <- stepp[nr]
                    }
                    
                  }
                  
                  d1 <- as.numeric(model.matrix[j, (NR1 + 1):(NR1 + 
                    NR2)])
                  model.matrix[j, (NR1 + 1):(NR1 + NR2)] <- sample(round(d1), 
                    NR2)
                }
                
                
                if (coin == 1) {
                  EDgenelist[index, 1] = j
                  EDgenelist[index, 2] = rownames(model.matrix)[j]
                  EDgenelist[index, 3] = fc
                  EDgenelist[index, 4] = log2(fc)
                } else {
                  EDgenelist[index, 1] = j
                  EDgenelist[index, 2] = rownames(model.matrix)[j]
                  EDgenelist[index, 3] = 1/fc
                  EDgenelist[index, 4] = log2(1/fc)
                }
                DEgeneID <- c(DEgeneID, rownames(model.matrix)[j])
                index <- index + 1
                
            } else {
                
                mean12 <- mean_fc_relation[j, 1]
                if (NDEfc == TRUE) {
                  fc <- mean_fc_relation[j, 2]
                } else {
                  fc = 1
                }
                
                if (fc >= 1) {
                  coin <- 0
                } else if (fc < 1) {
                  coin <- 1
                  fc <- 1/fc
                }
                
                {
                  mu = (coin * sqrt(fc) + (1 - coin)/sqrt(fc)) * 
                    mean12
                  tmp <- abs(mean_ref - mu)
                  idx <- order(tmp)[1]
                  
                  d1 <- as.numeric(data_norm[idx, ])
                  d1 <- round(d1 * mu/mean(d1))
                  
                  # model.matrix[j,1:NR1] <- sample(d1, NR1, replace=T)
                  stepp = as.numeric(quantile(d1, probs = seq(0, 
                    1, length.out = (NR1 + 1))))
                  
                  for (nr in 1:NR1) {
                    if (nr == 1) {
                      ind <- d1 >= stepp[nr] & d1 <= stepp[nr + 
                        1]
                    } else {
                      ind <- d1 > stepp[nr] & d1 <= stepp[nr + 
                        1]
                    }
                    
                    if (sum(ind) > 0) 
                      model.matrix[j, nr] <- sample(d1[ind], 
                        1, replace = TRUE) else {
                      model.matrix[j, nr] <- stepp[nr]
                    }
                  }
                  
                  d1 <- as.numeric(model.matrix[j, 1:NR1])
                  model.matrix[j, 1:NR1] <- sample(round(d1), 
                    NR1)
                  
                  
                }
                
                
                {
                  mu = ((1 - coin) * sqrt(fc) + coin/sqrt(fc)) * 
                    mean12
                  tmp <- abs(mean_ref - mu)
                  idx <- order(tmp)[1]
                  
                  d1 <- as.numeric(data_norm[idx, ])
                  d1 <- round(d1 * mu/mean(d1))
                  
                  # model.matrix[j,(NR1+1):(NR1+NR2)] <- sample(d1, NR2,
                  # replace=T)
                  stepp = as.numeric(quantile(d1, probs = seq(0, 
                    1, length.out = (NR2 + 1))))
                  
                  for (nr in 1:NR2) {
                    
                    if (nr == 1) {
                      ind <- d1 >= stepp[nr] & d1 <= stepp[nr + 
                        1]
                    } else {
                      ind <- d1 > stepp[nr] & d1 <= stepp[nr + 
                        1]
                    }
                    
                    if (sum(ind) > 0) 
                      model.matrix[j, nr + NR1] <- sample(d1[ind], 
                        1, replace = TRUE) else {
                      model.matrix[j, nr + NR1] <- stepp[nr]
                    }
                    
                  }
                  
                  d1 <- as.numeric(model.matrix[j, (NR1 + 1):(NR1 + 
                    NR2)])
                  model.matrix[j, (NR1 + 1):(NR1 + NR2)] <- sample(round(d1), 
                    NR2)
                  
                }
                
            }
            
            
        }
        
    }
    # Generate frequency matrix
    model.freq <- sweep(model.matrix, 2, apply(model.matrix, 
        2, sum), "/")
    
    if (model == "NegBinomial" || model == "ModelFree" || model == 
        "SingleCell") 
        model.combined <- round(sweep(model.freq, 2, libsizes, 
            "*"))
    
    if (model == "Multinom" || model == "Full") {
        # Generate multinomial distribution using frequen matrix
        model.combined <- matrix(0, nrow = EC, ncol = (NR1 + 
            NR2))
        rownames(model.combined) <- rownames(model.matrix)
        for (i in 1:(NR1 + NR2)) {
            model.combined[, i] <- rmultinom(n = 1, size = libsizes[i], 
                prob = model.freq[, i])
        }
    }
    return(list(count = model.combined, DiffAbundList = EDgenelist, 
        dataLabel = dataLabel))
}
