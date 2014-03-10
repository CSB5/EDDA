\name{computeAUC}
\alias{computeAUC}
\title{compute AUC values.}
\description{compute AUC values for each test.}
\usage{
computeAUC(obj,cutoff=1,numCores=10,
		DE.methods=c("Cuffdiff","DESeq","baySeq","edgeR","MetaStats","NOISeq"), 
		nor.methods=c("default","Mode","UQN","NDE")) 
}

\arguments{
  \item{obj}{
Object from testDATs().
}
  \item{cutoff}{
  cutoff for ROC curve. Default is 1.
}
  \item{numCores}{
Number of cores for parallelization. Default is 10.
}
  \item{DE.methods}{
Method list for differential abundance tests. Methods currently available include "Cuffdiff","DESeq", "baySeq","edgeR","MetaStats","NOISeq".
}
  \item{nor.methods}{
Normalization method list. Methods currently available include "default"(default normalization for each DE method), "Mode"(Mode normalization),"UQN"(Upper quartile normalization),"NDE"(non-differential expression normalization).
}
}

\references{
Luo Huaien, Li Juntao,Chia Kuan Hui Burton, Shyam Prabhakar, Paul Robson, 
Niranjan Nagarajan, The importance of study design for detecting 
differentially abundant features in high-throughput experiments, under review.
}
\author{Li Juntao, and Luo Huaien}
\examples{
data <- generateData(EntityCount=200)
test.obj <- testDATs(data,DE.methods="DESeq",nor.methods="default")
auc.obj  <- computeAUC(test.obj)
}
\keyword{ math }

