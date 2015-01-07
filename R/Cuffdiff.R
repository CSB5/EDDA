Cuffdiff <- function(counts, cond1, cond2){
  normFactors <- apply(counts, 2, sum) / 1e6;
  rpm <- sweep(counts, 2, normFactors, "/");	# normalize counts by per million read to get RPM
  
  results <- call_cuffdiff(rpm, normFactors, cond1, cond2);
  colnames(results) <- c("id", "RPM_A", "RPMvar_A", "RPM_B", "RPMvar_B", "pValue", "pAdj");	# setting colnames of results matrix

  return(results);
}


Cuffdiff_uqn <- function(counts, cond1, cond2){  
  returnList <- UQNnormalization(counts);
  uqn <- as.data.frame(returnList$normCounts); 
  normFactors <- returnList$normFactors;
  
  results <- call_cuffdiff(uqn, normFactors, cond1, cond2);
  colnames(results) <- c("id", "UQN_A", "UQNvar_A", "UQN_B", "UQNvar_B", "pValue", "pAdj");	# setting colnames of results matrix

  return(results);   
}


Cuffdiff_edda <- function(counts, conds, cond1, cond2, runID, winSize){
  returnList <- normalizeData(counts, conds, runID, winSize);
  normCounts <- as.data.frame(returnList$normCounts);
  normFactors <- returnList$normFactors;
  
  results <- call_cuffdiff(normCounts, normFactors, cond1, cond2);
  colnames(results) <- c("id", "counts_A", "var_A", "counts_B", "var_B", "pValue", "pAdj");	# setting colnames of results matrix

  return(results);
}


Cuffdiff_nde <- function(counts, cond1, cond2, runID, winSize){
  returnList <- normalizeNDE(counts, runID, winSize);
  normCounts <- as.data.frame(returnList$normCounts);
  normFactors <- returnList$normFactors;
  
  results <- call_cuffdiff(normCounts, normFactors, cond1, cond2);
  colnames(results) <- c("id", "counts_A", "var_A", "counts_B", "var_B", "pValue", "pAdj");	# setting colnames of results matrix

  return(results);
}


call_cuffdiff <- function(counts, normFactors, cond1, cond2){
	results <- data.frame(matrix(nrow=length(counts[,1]), ncol=7));	# create empty data frame to hold results
	for (i in seq(from=1, to=length(counts[,1]), by=1)){
		id <- rownames(counts)[i];
		#rownames(results)[i] <- rownames(counts)[i];
		
		if(cond1 == 1){
			mean_A <- counts[i,1:cond1];
		} else{
			mean_A <- rowMeans(counts[i,1:cond1]);
		}
		
		if(cond2 == 1){
			mean_B <- counts[i,(cond1+1):(cond1+cond2)];
		} else{
			mean_B <- rowMeans(counts[i,(cond1+1):(cond1+cond2)]);
		}
		
		if(cond1 == 1 || cond2 == 1){
			if(cond1 == 1){
				var_A = var_B = var(as.integer(counts[i,(cond1+1):(cond1+cond2)]));
			} else{
				var_A = var_B = var(as.integer(counts[i,]));	
			}			
		} else{
			var_A <- var(as.integer(counts[i,1:cond1]));
			var_B <- var(as.integer(counts[i,(cond1+1):(cond1+cond2)]));
		}
		
		pVal <- .Call("cuffdiff_wrapper", mean_A, var_A, mean_B, var_B)
		results[i,1:6] <- c(id, mean_A, var_A, mean_B, var_B, pVal);	# storing results in table
	} # end for loop
	results[,7] <- p.adjust(results[,6], method = "BH");
	
	return(results);
}

