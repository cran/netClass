\name{cv.pac}
\alias{cv.pac}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cross validation for Pathway Activities Classification(PAC)
}
\description{
Cross validation for Pathway Activities Classification(PAC) using Logistic regression 
model for classification. Implementation of the Pathway Activities Classification  
by CROG algorithm.
}
\usage{
cv.pac(x=x, y=y, folds=10, repeats=5, parallel = TRUE, cores = NULL, 
		DEBUG=TRUE, Gsub=matrix(1,100,100), seed=1234)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
 a p x n matrix of expression measurements with p samples and n genes.
}
  \item{y}{
a factor of length p comprising the class labels.
}
  \item{folds}{
number of -folds cross validation (CV)
}
  \item{repeats}{
number of CV repeat times
}
  \item{parallel}{
paralle computing or not
}
  \item{cores}{
cores used in parallel computing
}
  \item{DEBUG}{
show debugging information in screen or not.
}
 
  \item{Gsub}{
Adjacency matrix of Protein-protein intersction network
}
  \item{seed}{
seed for random sampling.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
a LIST for Cross-Validation results
\item{auc }{The AUC values of each test fold}
\item{fits }{The tranined models for traning folds}
\item{feat }{The feature selected by each by the fits }
\item{labels }{the original lables for training}

}
\references{
Lee E, Chuang H-Y, Kim J-W, Ideker T, Lee D (2008) Inferring Pathway Activity toward Precise Disease Classification. PLoS Comput Biol 4(11): e1000217. 
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
\examples{
library(netClass)
 
data(expr)
data(ad.matrix)
x <- expr$genes
y <- expr$y

library(KEGG.db)
r.pac <- cv.pac(x=x, y=y, folds=3, repeats=1, parallel=FALSE, cores=2, DEBUG=TRUE,
				Gsub=ad.matrix,seed=1234)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ pac }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
