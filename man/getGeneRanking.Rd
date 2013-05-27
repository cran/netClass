\name{getGeneRanking}
\alias{getGeneRanking}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get gene ranking based on geneRank algorithm.
}
\description{
Get the ranking of differential expression of genes on graph using geneRank algorithm.
}
\usage{
getGeneRanking(x = x, y = y, Gsub = Gsub, d = d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
gene expression data
}
  \item{y}{
class labels
}
  \item{Gsub}{
Adjacency matrix of Protein-protein intersction network
}
  \item{d}{
damping factor for GeneRank, defaults value is 0.5
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{r }{ranking of each gebes on graph}
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
See Also as {pGeneRank}
}
\examples{
library(netClass)
data(expr)
data(ad.matrix)
ex.sum <- expr$genes
y <- expr$y

#r= getGeneRanking(x = ex.sum, y = y, Gsub = ad.matrix, d = 0.5)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ FrSVM  }
%\keyword{ GeneRank }% _
