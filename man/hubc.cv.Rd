\name{hubc.cv}
\alias{hubc.cv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cross validation for hub nodes classification
}
\description{
Cross validation for hub nodes classification, which described in Taylor et al.(2009).

}
\usage{
hubc.cv(x, y, folds=10, repeats=5, parallel = TRUE, cores = NULL, DEBUG=TRUE, nperm=500,
		node.ct=0.98, Gsub=matrix(1,100,100), Gs=Gs,seed=1234,Cs=10^c(-3:3))
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
  \item{nperm}{
number of permutation test steps
}
  \item{node.ct}{
cut off value for select highly quantile nodes in a nwtwork. Defaults to \code{0.98)}.
}
  \item{Gsub}{
an adjacency matrix that represents the underlying biological network.
}
  \item{Gs}{
Undirected of graph with adjacency matrix Gsub.
}
  \item{seed}{
Seed for random sampling.
}
  \item{Cs}{
Soft-margin tuning parameter of the SVM. Defaults to \code{10^c(-3:3)}.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{auc }{The AUC values of each test fold}
\item{fits }{The tranined models for traning folds}
\item{feat }{The selected features of each training folds}
\item{labels }{the original lables for training}
%% ...
}
\references{
Taylor et al.(2009)Dynamic modularity in protein interaction networks predicts breast cancer outcome, Nat. Biotech.: doi: 10.1038/nbt.1522
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
 data(Gs2)
 data(ad.matrix)
 x <- expr$genes
 y <- expr$y
r.hubC <- hubc.cv(x=x, y=y, folds=3, repeats=1, parallel=FALSE, cores=2, DEBUG=TRUE,
				nperm=2, Gsub=ad.matrix,Gs=Gs2,node.ct=0.6,Cs=10^(-3:3))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hubc }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
