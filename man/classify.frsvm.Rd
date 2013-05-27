\name{classify.frsvm}
\alias{classify.frsvm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Training and predicting using FrSVM
}
\description{
Training and predicting using FrSVM
}
\usage{
classify.frsvm(fold, cuts, x, y, cv.repeat, DEBUG = DEBUG, Gsub = Gsub, d = d, op = op,
				aa = aa, Cs = Cs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fold}{
number of folds to perform
}
  \item{cuts}{
list for randomly divide the training set in to x-x-CV
}
   \item{x}{
expression data 
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
  \item{d}{
damping factor for GeneRank, defaults value is 0.5
}
  \item{op}{
 the uper bound of top ranked genes
}
  \item{aa}{
the lower bound of top ranked genes
}
  \item{Cs}{
soft-margin tuning parameter of the SVM. Defaults to \code{10^c(-3:3)}.
}
}
%\details{
%%  ~~ If necessary, more details than the description in \name{FrSVM.cv}
%}
\value{
%%  ~Describe the value returned
\item{fold }{the recored for test fold}
\item{auc }{The AUC values of test fold}
\item{train }{The tranined models for traning folds}
\item{feat }{The feature selected by each by the train}
%% ...
}
\references{
Yupeng Cun, Holger Frohlich (2012) Integrating Prior Knowledge Into Prognostic Biomarker Discovery Based on Network Structure.arXiv:1212.3214 \cr
Winter C, Kristiansen G, Kersting S, Roy J, Aust D, et al. (2012) Google Goes Cancer: Improving Outcome Prediction for Cancer Patients by Network-Based Ranking of Marker Genes. PLoS Comput Biol 8(5): e1002511. doi:10.1371/journal.pcbi.1002511 \cr
}
\author{
Yupeng Cun \email{yupeng.cun@gmail.com}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See Also as {FrSVM.cv}
}
\examples{
#see \name{FrSVM.cv}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ FrSVM }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
