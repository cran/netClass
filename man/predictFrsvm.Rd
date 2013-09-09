\name{predictFrsvm}
\alias{predictFrsvm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predicting the test data using frsvm trained model
}
\description{
Predicting the test data using frsvm trained model
}
\usage{
predictFrsvm(train = train, x = x, y = y, DEBUG = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{train}{
trained model
}
  \item{x}{
expression data for testing
}
  \item{y}{
  class labels
}
  \item{DEBUG}{
show debugging information in screen more or less.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
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
See Also as cv.frsvm
}
\examples{
#see cv.frsvm
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ frsvm }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
