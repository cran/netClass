\name{train.hubc}
\alias{train.hubc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Predicting the data using hub nodes classification model
}
\description{
 Predicting the data using hub nodes classification model
}
\usage{
train.hubc(x = x, y = y, DEBUG = FALSE, Gsub = Gsub, gHub = gHub, 
		hubs = hubs, nperm = 500, node.ct = 0.95, Cs = 10^(-3:3))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
gene expression data for training.
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
The list returned
%%  If it is a LIST, use
\item{trained }{The tranined models for traning folds}
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
See Also as cv.hubc
}
\examples{
#See cv.hubc
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hubc }
