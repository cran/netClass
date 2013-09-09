\name{predictPac}
\alias{predictPac}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predicting the test data using pac trained model
}
\description{
Predicting the test data using pac trained model
}
\usage{
predictPac(train = train, x = x, y = y, int = int, DEBUG = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{train}{
%%     ~~Describe \code{train} here~~
}
  \item{x}{
gene expression data for the testing data
}
  \item{y}{
a factor of length p comprising the class labels.
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
See Also as cv.pac
}
\examples{
#see cv.pac
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ pac }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
