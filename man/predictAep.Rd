\name{predictAep}
\alias{predictAep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predicting the test tdata using aep trained model
}

\description{
Predicting the test data using aep trained model
}
\usage{
predictAep(train = train, x, y, DEBUG = FALSE, Gsub = Gsub)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{train}{
trained model
}
  \item{x}{
gene expression data for testing
}
  \item{y}{
class labels
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
The value returned
%%  If it is a LIST, use
\item{auc }{The AUC values of test fold}
%%  \item{comp2 }{Description of 'comp2'}
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
See Also as cv.aep
}
\examples{
#see cv.aep
}
\keyword{ aep }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
