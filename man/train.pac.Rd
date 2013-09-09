\name{train.pac}
\alias{train.pac}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Training the data using pac methods
}
\description{
Training the data using pac methods
}
\usage{
train.pac(x = x, y = y, int = int, DEBUG = FALSE, Gsub = Gsub)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
gene expression data for the training data
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
  \item{Gsub}{
 an adjacency matrix that represents the underlying biological network.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
 the value returned
%%  If it is a LIST, use
\item{trained }{The tranined models for traning folds}
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
