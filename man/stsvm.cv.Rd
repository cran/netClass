\name{stsvm.cv}
\alias{stsvm.cv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cross validation for smoothed t-statistic to select significant top ranked differential expressed genes
}
\description{
Cross validation for smoothed t-statistic to select significant top ranked differential expressed genes
}
\usage{
stsvm.cv(x=x, x.mi=NULL,y=y, folds=5,Gsub=matrix(1,100,100),repeats=3, parallel=FALSE,
	cores=2, DEBUG=TRUE, pt.pvalue = 0.05, op = 0.85, pt.step = 1000, a = 1, p=2,
	allF=TRUE, seed=1234, Cs = 10^c(-3:3) )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A p x n matrix of expression measurements with p samples and n genes.
}
  \item{x.mi}{
A p x m matrix of expression measurements with p samples and m miRNAs.
}
  \item{y}{
A factor of length p comprising the class labels.
}
  \item{folds}{
Folds number of folds to perform
}
  \item{Gsub}{
An adjacency matrix that represents the underlying biological network.
}
  \item{repeats}{
Number of how often to repeat the x-fold cross-validation
}
  \item{parallel}{
Use parallel computing or not
}
  \item{cores}{
Number of cores will used when parallel is TRUE
}
  \item{DEBUG}{
Show debugging information in screen more or less.
}
  \item{pt.pvalue}{
Cut off p-value of permutation test
}
  \item{op}{
Pptimal on top op% ranked genes
}
  \item{pt.step}{
permutation test steps
}
  \item{a}{
constant value of random walk kernel
}
  \item{p}{
random walk step(s) of random walk kernel
}
  \item{allF}{
Using all features (TRUE) or only these genes mapped to prior information (FALSE). 
}
  \item{seed}{
seed for random sampling.
}
  \item{Cs}{
Soft-margin tuning parameter of the SVM. Defaults to \code{10^c(-3:3)}.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
a LIST for Cross-Validation results
\item{auc }{The AUC values of each test fold}
\item{fits }{The tranined models for traning folds}
\item{feat }{The feature selected by each by the fits }
\item{labels }{the original lables for training}
%% ...auc=auc, fits=fits, feat= genes, labels=y
}

\references{
Yupeng Cun, Holger Frohlich (2013) Network and Data Integration for Biomarker Signature Disvovery via Network Smoothed T-Statistics
}
\author{
Yupeng Cun \email{yupeng.cun@gmail.com}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
library(netClass)
data(expr)
data(ad.matrix)
x <- expr$genes
y <- expr$y

r.stsvm <- stsvm.cv(x=x[,1:100],x.mi=NULL,y=y,folds=3,Gsub=ad.matrix,repeats=1, parallel=FALSE,
			 cores=2,DEBUG=TRUE,pt.pvalue=0.05,op=0.9,pt.step=3,a=1,p=2,allF=TRUE,
			 seed=1234,Cs=10^(-3:3))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ stSVM }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
