\name{train.frsvm}
\alias{train.frsvm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Training the data using frsvm method
}
\description{
Training the data using frsvm methods
}
\usage{
train.frsvm(x = x, y = y, DEBUG = FALSE, Gsub = Gsub, d = 0.85, op
                 = 10, aa = 50, Cs = 10^(-3:3))
          }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Expression data for training
}
  \item{y}{
Class labels
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
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
The value list returned
%%  If it is a LIST, use
\item{train }{The tranined models for traning folds}
\item{feat }{The feature selected by each by the train}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Yupeng Cun \email{yupeng.cun@gmail.com}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See Also as cv.frsvm
}
\examples{
#see cv.frsvm
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ frsvm }
\keyword{ GeneRank }% __ONLY ONE__ keyword per line
