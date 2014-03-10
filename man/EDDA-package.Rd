\name{EDDA-package}
\alias{EDDA-package}
\alias{EDDA}
\docType{package}
\title{
Experimental Design in Differential Abundance analysis
}
\description{
EDDA aids in the design of a range of common experiments including RNA-seq, 
Nanostring assays, RIP-seq and Metagenomic sequencing, and enables researchers 
to comprehensively investigate the impact of experimental decisions on the 
ability to detect differential abundance.
}
\details{
\tabular{ll}{
Package: \tab EDDA\cr
Type: \tab Package\cr
Version: \tab 0.99.2\cr
Date: \tab 2014-02-12\cr
License: \tab GPL (>= 2)\cr
}
generateData()
testDATs()
computeAUC()
plotROC()
plotPRC()
}
\author{Li Juntao, Luo Huaien, Chia Kuan Hui Burton, Niranjan Nagarajan

Maintainer: Li Juntao<lij9@gis.a-star.edu.sg>, 
            Luo Huaien<luoh2@gis.a-star.edu.sg>, 
	    Niranjan Nagarajan <nagarajann@gis.a-star.edu.sg>
}
\references{
Luo Huaien, Li Juntao,Chia Kuan Hui Burton, Shyam Prabhakar, Paul Robson, 
Niranjan Nagarajan, The importance of study design for detecting 
differentially abundant features in high-throughput experiments, under review.
}
\keyword{ package }

\examples{
data <- generateData(EntityCount=500)
test.obj <- testDATs(data,DE.methods=c("DESeq","edgeR"),nor.methods="default")
auc.obj  <- computeAUC(test.obj)
plotROC(auc.obj)
plotPRC(auc.obj)
}
