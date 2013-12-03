\name{cv.aep}
\alias{cv.aep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cross validation for aepSVM (aepSVM)
}
\description{
Cross validation for aepSVM (aepSVM) using SAM to select significant differential expressed genes
}
\usage{
cv.aep(x, y, folds = 10, repeats = 5, parallel = FALSE, cores
                 = 2, DEBUG = TRUE, Gsub = matrix(1, 100, 100), 
                 Cs = 10^(-3:3), seed = 1234)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a p x n matrix of expression measurements with p samples and n genes.
}
  \item{y}{
a factor of length p comprising the class labels.
}
  \item{folds}{
number of -folds cross validation (CV)
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
 
  \item{Gsub}{
Adjacency matrix of Protein-protein interaction network
}
  \item{Cs}{
soft-margin tuning parameter of the SVM. Defaults to \code{10^c(-3:3)}.
}
  \item{seed}{
seed for random sampling.
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
\item{feat }{The feature selected by each by the fits}
\item{labels }{the original lables for training}
%% ...
}
\references{
Guo et al., Towards precise classification of cancers based on robust gene functional expression profiles. BMC Bioinformatics 2005, 6:58. 
}
\author{
 Yupeng Cun \email{yupeng.cun@gmail.com}
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% See Also as \code{cv.aep}
%}
\examples{
library(netClass)
data(expr)
data(ad.matrix)
x <- expr$genes
y <- expr$y

 library(KEGG.db)
#r.aep <- cv.aep(x[,1:500], y, folds=3, repeats=1, parallel=FALSE,cores=2,
#			Gsub=ad.matrix,	Cs=10^(-3:3),seed=1234,DEBUG=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ aep }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
