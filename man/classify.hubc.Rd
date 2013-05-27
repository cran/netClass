\name{classify.hubc}
\alias{classify.hubc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Training and predicting using hub nodes classification methods
}
\description{
 Training and predicting using hub nodes classification methods
}
\usage{
classify.hubc(fold, r, cuts, x, y, cv.repeat, DEBUG = DEBUG, gHub = gHub, hubs = hubs,
				nperm = nperm, node.ct = node.ct, Cs = Cs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fold}{
number of -fold cross validation (CV)
}

  \item{cuts}{
list for randomly divide the training set in to x-x-fold CV
}

  \item{x}{
gene expression data.
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
 
  \item{r}{
repeat order for CV
}
  \item{gHub}{
Subgraph of hubs of graph Gs
}
  \item{hubs}{
Hubs in graph Gs
}
  \item{nperm}{
number of permutation test steps
}
  \item{node.ct}{
cut off value for select highly quantile nodes in a nwtwork. Defaults to \code{0.98)}.
}
  \item{Cs}{
Soft-margin tuning parameter of the SVM. Defaults to \code{10^c(-3:3)}.
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
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hubc }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
