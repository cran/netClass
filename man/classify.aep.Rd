\name{classify.aep}
\alias{classify.aep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Training and predicting using aepSVM (aepSVM) classification methods
}
\description{
Training and predicting using aepSVM (aepSVM) classification methods
}
\usage{
classify.aep(fold, cuts, Cs, x, y, cv.repeat, int, DEBUG = DEBUG, Gsub)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fold}{
number of -folds cross validation (CV)
}
  \item{cuts}{
list for randomly divide the training set in to x-x-folds CV
}

  \item{Cs}{
soft-margin tuning parameter of the SVM. Defaults to \code{10^c(-3:3)}.
}
  \item{x}{
gene expression data
}
  \item{y}{
class labels
}

  \item{cv.repeat}{
model for one CV training and predicting
}
  \item{int}{
Intersect of genes in network and gene expression profile. 
}
  \item{DEBUG}{
show debugging information in screen more or less.
}

  \item{Gsub}{
 an adjacency matrix that represents the underlying biological network.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{fold }{the recored for test fold}
\item{auc }{The AUC values of test fold}
\item{train }{The tranined models for traning folds}
\item{feat }{The feature selected by each by the train}
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

\seealso{
See Also as \code{aep.cv}
}
\examples{
#see \code{aep.cv}


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ aep }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
