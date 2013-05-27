\name{classify.pac}
\alias{classify.pac}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Training and predicting using PAC classification methods
}
\description{
Training and predicting using PAC classification methods
}
\usage{
classify.pac(fold, cuts, x, y, cv.repeat, int, DEBUG = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fold}{
number of -folds cross validation (CV)
}
  \item{cuts}{
list for randomly divide the training set in to x-x-folds CV
}

  \item{x}{
gene expression data
}
  \item{y}{
a factor of length p comprising the class labels.
}
  \item{cv.repeat}{
model for one CV training and predicting
}

 \item{int}{
Intersect of genes in network and gene expression profile. 
}
  \item{DEBUG}{
show debugging information in screen or not.
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
Lee E, Chuang H-Y, Kim J-W, Ideker T, Lee D (2008) Inferring Pathway Activity toward Precise Disease Classification. PLoS Comput Biol 4(11): e1000217. doi:10.1371/journal.pcbi.1000217
}
\author{
Yupeng Cun \email{yupeng.cun@gmail.com}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See Also as {pac.cv}
}
\examples{
#see \name{FrSVM.cv}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ pac }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
