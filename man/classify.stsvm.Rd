\name{classify.stsvm}
\alias{classify.stsvm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Training and predicting using stSVM classification methods
}
\description{
Training and predicting using stSVM classification methods
}
\usage{
classify.stsvm(fold, cuts, ex.sum, x, p, a, y, cv.repeat, DEBUG = DEBUG, Gsub=Gsub,
		op = op, aa = aa, dk = dk, dk.tf = dk.tf, seed = seed, Cs = Cs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fold}{
number of folds to perform
}
  \item{cuts}{
list for randomly divide the training set in to x-x-folds CV
}

  \item{ex.sum}{
expression data 
}
  \item{x}{
expression data 
}
  \item{a}{
constant value of random walk kernel
}
  \item{p}{
random walk step(s) of random walk kernel
}
  \item{y}{
a factor of length p comprising the class labels.
}
  \item{cv.repeat}{
model for one CV training and predicting
}
  \item{DEBUG}{
show debugging information in screen more or less.
}
  \item{Gsub}{
 an adjacency matrix that represents the underlying biological network.
}
  \item{op}{
optimal on top op% ranked genes
}
  \item{aa}{
permutation test steps
}
  \item{dk}{
Random Walk Kernel matrix of network
}
  \item{dk.tf}{
cut off p-value of permutation test
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
%%  If it is a LIST, use
\item{fold }{the recored for test fold}
\item{auc }{The AUC values of test fold}
\item{train }{The tranined models for traning folds}
\item{feat }{The feature selected by each by the train }
%% ...
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
#see \name{stsvm.cv}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ stSVM }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line%
