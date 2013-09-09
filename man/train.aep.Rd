\name{train.aep}
\alias{train.aep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Training the data using aep methods
}
\description{
Training the data using aep methods
}
\usage{
train.aep(x = x, y = y, DEBUG = FALSE, int = int, Gsub = Gsub, Cs = 10^(-3:3))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
expression data for training
}
  \item{y}{
a factor of length p comprising the class labels.
}
  \item{DEBUG}{
show debugging information in screen more or less.
}
  \item{int}{
Intersect of genes in network and gene expression profile. 
}
  \item{Gsub}{
 an adjacency matrix that represents the underlying biological network.
}
  \item{Cs}{
soft-margin tuning parameter of the SVM. Defaults to \code{10^c(-3:3)}.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
 The returned lists
\item{trained }{The tranined models for traning folds}
\item{sig.genes }{The differential expressed feature}
%% ...
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
See Also as cv.aep
}
\examples{
#see cv.aep

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ aep }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
