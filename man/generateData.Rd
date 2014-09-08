\name{generateData}
\alias{generateData}
\title{generate count data}
\description{Simulate count data using different models and settings.}
\usage{
generateData(SimulModel="Full", SampleVar="medium", 
		ControlRep=5, CaseRep=ControlRep, EntityCount=1000, FC="Norm(2,1)",
		perDiffAbund=0.1, upPDA=perDiffAbund/2, downPDA=perDiffAbund/2,
		numDataPoints=100, AbundProfile = "HBR", modelFile = NULL, minAbund=10,varLibsizes=0.1,
		inputCount=NULL,inputLabel=NULL,SimulType="auto")
}
\arguments{
  \item{SimulModel}{
Simulation model used. Default is "Full".
SimulModel="NegBinomial" is negtive binomial model which generates data using negtive binomial distribution.
SimulModel="Multinom" is multinomial model which generates data mimicking the multinomial sampling process.
SimulModel="Full" is a model which combines "NegBinomial" and "Multinom".
SimulModel="ModelFree" uses model free approach to generate data by sub-sampling counts from modelFile (if modelFile != NULL) or from input File (if modelFile == NULL)  
}
  \item{SampleVar}{
Sample variation: Default is "medium".
It could be "low", "medium" and "high" or a real number. 
}
  \item{ControlRep}{
Number of replicates for control group. Default is 5.
}
  \item{CaseRep}{
Number of replicates for case group. Default is same as ControlRep.
}
  \item{EntityCount}{
Entity count.  Default is 1000.
}
  \item{FC}{
Fold change type. It can be "Norm(mu,sigma)", "logNorm(mu,sigma)", "log2Norm(mu,sigma)" or "Unif(a,b)". mu,sigma and a,b need be predefined. Default is "Norm(2,1)".
}
  \item{perDiffAbund}{
Percentage of differential abundance. Default is 0.1 (i.e "10 percent").
}
  \item{upPDA}{
Percentage of up-regulated differential abundance. Default is perDiffAbund/2.
}
  \item{downPDA}{
Percentage of down-regulated differential abundance. Default is perDiffAbund/2.
}
  \item{numDataPoints}{
Number of data points. Default is 100.
}
  \item{AbundProfile}{
AbundProfile for average abundance profile. It can be either the different profiles used in the paper ("HBR", "BP" and "Wu") or it can be location of the abundance profile. Default is "HBR".
}

 \item{modelFile}{
Sample data file for model free approach. Default is NULL.
If modelFile = NULL, Model Free approach will subsample from Input file.
If modelFile = "SingleCell", Model Free approach will subsample from the available single cell RNA-seq data.
if modelFile is the name of a count file, this count file will be used as sample file for sub-sampling.
}

  \item{minAbund}{
Minimum abundance cutoff. Default is 10.
}
  \item{varLibsizes}{
Variability between library sizes. Default is 0.1.
}

\item{outlier}{
Percentage of added outliers. Default is 0.0 (i.e. no outlier).
}
  \item{inputCount}{
Input count file. Default is NULL. If not NULL, it learns the parameters (modelFile, SampleVar, perDiffAbund, upPDA, and downPDA) from count data.
}
  \item{inputLabel}{
Label of input count file. The label should be sequence of 0 or 1. Default is NULL.
}
  \item{SimulType}{
Simulation type. It is used only when user has pilot data. Default is "auto".
SimulType = "auto", all the parameters (EntityCount, ControlRep, CaseRep, numDataPoint and others) are learned from user's pilot data. 
SimulType = "auto1", ControlRep, CaseRep, numDataPoint are specified by user input; while EntityCount and all others are learned from user's pilot data.
SimulType = "auto2", EntityCount, ControlRep, CaseRep, numDataPoint are specified by user input; while all others are learned from user's pilot data.
}

}

\value{
\item{count }{Count matrix.}
\item{DiffAbundList }{Differential abundance list.}
\item{dataLabel }{Data label.}
}
\references{
Luo Huaien, Li Juntao,Chia Kuan Hui Burton, Shyam Prabhakar, Paul Robson, 
Niranjan Nagarajan, The importance of study design for detecting 
differentially abundant features in high-throughput experiments, under review.
}
\author{Li Juntao, Luo Huaien, Chia Kuan Hui Burton, Niranjan Nagarajan}
\examples{
# generate data with all default options.
data <- generateData()
dim(data$count)
dim(data$DiffAbundList)
data$dataLabel

# generate data with input count.
x <- matrix(rnbinom(1000*15,size=1,mu=10), nrow=1000, ncol=15);
x.lable=c(rep(0,10),rep(1,5))
x[1:50,11:15] <- x[1:50,11:15]*10
x.name=paste("g",1:1000,sep="");
write.table(cbind(x.name,x),"count.txt",row.names =FALSE, sep ='\t')

data <- generateData(inputCount="count.txt",inputLabel=x.lable)
dim(data$count)
dim(data$DiffAbundList)
data$dataLabel

# or generate data with input count and redefined parameters.
data <- generateData(inputCount="count.txt",inputLabel=x.lable,
                     ControlRep=10,CaseRep=10,EntityCount=3000,SimulType="auto2")
dim(data$count)
dim(data$DiffAbundList)
data$dataLabel

}
\keyword{design}

