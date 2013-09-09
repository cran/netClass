\name{cv.frsvm}
\alias{cv.frsvm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cross validation for FrSVM
}
\description{
Cross validation for FrSVM, an R algorithm, which integrates protein-protein interaction network information into gene selection for microarry classification
}
\usage{
cv.frsvm(x, y, folds = 10, Gsub = matrix(1, 100, 100), repeats
                 = 5, parallel = FALSE, cores = 2, DEBUG = FALSE, d =
                 0.85, top.uper = 10, top.lower = 50, seed = 1234, Cs =
                 10^c(-3:3))

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
gene expression data
}
  \item{y}{
class labels
}
  \item{folds}{
number of -folds cross validation (CV)
}
  \item{Gsub}{
Adjacency matrix of Protein-protein intersction network
}
  \item{repeats}{
number of CV repeat times
}
  \item{parallel}{
paralle computing or not
}
  \item{cores}{
cores used in parallel computing
}
  \item{DEBUG}{
show more results or not
}
  \item{d}{
damping factor for GeneRank, defaults value is 0.5
}
  \item{top.uper}{
the uper bound of top ranked genes
}
  \item{top.lower}{
the lower bound of top ranked genes
}
  \item{seed}{
Seed for random sampling.
}
  \item{Cs}{
soft-margin tuning parameter of the SVM. Defaults to \code{10^c(-3:3)}.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
a LIST for Cross-Validation results
\item{auc }{The AUC values of each test fold}
\item{fits }{The tranined models for traning folds}
\item{feat }{The feature selected by each by the fits}
\item{labels }{the original lables for training}

}
\references{
Yupeng Cun, Holger Frohlich (2012) Integrating Prior Knowledge Into Prognostic 
	Biomarker Discovery Based on Network Structure, arXiv:1212.3214 \cr 
Winter C, Kristiansen G, Kersting S, Roy J, Aust D, et al. (2012) Google Goes 
	Cancer: Improving Outcome Prediction for Cancer Patients by Network-Based 
	Ranking of Marker Genes. PLoS Comput Biol 8(5): e1002511. \cr

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
### 
r.frsvm <-cv.frsvm(x[,1:200], y, folds=3,Gsub=ad.matrix, repeats=1, parallel=FALSE,
	cores=2, DEBUG=TRUE,d=.85,top.uper=5,top.lower=15,seed=1234,Cs=10^c(-3:3))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ frsvm}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
