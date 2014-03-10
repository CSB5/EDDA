\name{testDATs}
\alias{testDATs}
\title{Run differential abundance testings}
\description{Perform differential abundance testing on simulated count data.}
\usage{
testDATs(data, numCores=10, minCountsThreshold=0,
	DE.methods=c("Cuffdiff","DESeq","baySeq","edgeR","MetaStats","NOISeq"), 
	nor.methods=c("default","Mode","UQN","NDE"),method.list=NULL)
}

\arguments{
  \item{data}{ Data object from generateData() function or predifined data 
  object similar to the output of generateData().
}

  \item{numCores}{
Number of cores for parallelization. Default is 10.
}
  \item{minCountsThreshold}{
Minimum counts threshold for filtering. Default is 0.
}
  \item{DE.methods}{
Method list for differential expression tests. Methods currently available 
include  "Cuffdiff","DESeq","baySeq","edgeR","MetaStats","NOISeq".
}
  \item{nor.methods}{
Normalization method list. Methods currently available include "default"
(default normalization for each DE method),"Mode"(Mode normalization),
"UQN"(Upper quartile normalization),
"NDE"(non-differential expression normalization).
}
\item{method.list}{
The method list for the combination of DE.methods and nor.methods. 
Default is NULL.
}

}

\value{
\item{data }{Data object from generateData() function.}
\item{filterCounts }{filtered count data.}
\item{Cuffdiff}{Result form Cuffdiff with default normalization.}
\item{Cuffdiff_uqn}{Result form Cuffdiff with Upper quartile normalization normalization.}
\item{Cuffdiff_Mode}{Result form Cuffdiff with Mode normalization.}
\item{Cuffdiff_nde}{Result form Cuffdiff with non-differential expression normalization.}
\item{DESeq}{Result form DESeq with default normalization.}
\item{DESeq_uqn}{Result form DESeq with Upper quartile normalization normalization.}
\item{DESeq_Mode}{Result form DESeq with Mode normalization.}
\item{DESeq_nde}{Result form DESeq with non-differential expression normalization.}
\item{baySeq}{Result form baySeq with default normalization.}
\item{baySeq_uqn}{Result form baySeq with Upper quartile normalization normalization.}
\item{baySeq_Mode}{Result form baySeq with Mode normalization.}
\item{baySeq_nde}{Result form baySeq with non-differential expression normalization.}
\item{edgeR}{Result form edgeR with default normalization.}
\item{edgeR_uqn}{Result form edgeR with Upper quartile normalization normalization.}
\item{edgeR_Mode}{Result form edgeR with Mode normalization.}
\item{edgeR_nde}{Result form edgeR with non-differential expression normalization.}
\item{MetaStats}{Result form MetaStats with default normalization.}
\item{MetaStats_uqn}{Result form MetaStats with Upper quartile normalization normalization.}
\item{MetaStats_Mode}{Result form MetaStats with Mode normalization.}
\item{MetaStats_nde}{Result form MetaStats with non-differential expression normalization.}
\item{NOISeq}{Result form NOISeq with default normalization.}
\item{NOISeq_uqn}{Result form NOISeq with Upper quartile normalization normalization.}
\item{NOISeq_Mode}{Result form NOISeq with Mode normalization.}
\item{NOISeq_nde}{Result form NOISeq with non-differential expression normalization.}
}

\references{
Luo Huaien, Li Juntao,Chia Kuan Hui Burton, Shyam Prabhakar, Paul Robson, 
Niranjan Nagarajan, The importance of study design for detecting 
differentially abundant features in high-throughput experiments, under review.
}
\author{Li Juntao, Luo Huaien, Chia Kuan Hui Burton, Niranjan Nagarajan}

\examples{

data <- generateData(EntityCount=100)
test.obj <- testDATs(data,nor.methods="default")
test.obj <- testDATs(data,DE.methods="DESeq")


# test data with input count.
x <- matrix(rnbinom(1000*15,size=1,mu=10), nrow=1000, ncol=15);
x[1:50,11:15] <- x[1:50,11:15]*10
x.name=paste("g",1:1000,sep="");
write.table(cbind(x.name,x),"count.txt",row.names =FALSE, sep ='\t')

x <- read.table("count.txt",head=TRUE,sep='\t')
x.count <- x[,2:16]
x.lable=c(rep(0,10),rep(1,5))
row.names(x.count) <- x[,1]
data <- list(count=x.count,dataLabel=x.lable)
test.obj <- testDATs(data,DE.methods=c("DESeq","edgeR"),nor.methods="default")

}
\keyword{methods }
