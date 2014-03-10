\name{plotPRC}
\alias{plotPRC}
\title{plot precision-recall curves}
\description{plot precision-recall curves for each test.}
\usage{
plotPRC(obj,DE.methods=c("Cuffdiff","DESeq","baySeq","edgeR","MetaStats","NOISeq"), 
		nor.methods=c("default","Mode","UQN","NDE"),
		plot_type = "o",plot_pch = 20,plot_lwd = 1.75,plot_cex = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
Object from testDATs().
}
  \item{DE.methods}{
Method list for differential expression tests. Methods currently available include "Cuffdiff","DESeq","baySeq","edgeR","MetaStats","NOISeq".
}
  \item{nor.methods}{
Normalization method list. Methods currently available include "default"(default normalization for each DE method),"Mode"(Mode normalization),"UQN"(Upper quartile normalization),"NDE"(non-differential expression normalization).
}
  \item{plot_type}{type option in plot.}
  \item{plot_pch}{pch option in plot.}
  \item{plot_lwd}{lwd option in plot.}
  \item{plot_cex}{cex option in plot.}

}

\references{
Luo Huaien, Li Juntao,Chia Kuan Hui Burton, Shyam Prabhakar, Paul Robson, 
Niranjan Nagarajan, The importance of study design for detecting 
differentially abundant features in high-throughput experiments, under review.
}
\author{Li Juntao and Chia Kuan Hui Burton}
\examples{
data <- generateData(EntityCount=500)
test.obj <- testDATs(data,DE.methods=c("DESeq","edgeR"),nor.methods="default")
auc.obj  <- computeAUC(test.obj)
plotPRC(auc.obj)
}
\keyword{dplot}
