\name{pGeneRANK}
\alias{pGeneRANK}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GeneRANK
}
\description{
Ranking gene based on Googles's PageRank algorithm
}
\usage{
pGeneRANK(W, ex, d, max.degree = Inf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{W}{
adjacency matrix of graph
}
  \item{ex}{
the fold change/ diffiencicial expression of genes
}
  \item{d}{
damping factor for GeneRank, defaults value is 0.5
}
  \item{max.degree}{
Max degree of graph 
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{r }{ranking of each gebes on graph}

%% ...
}
\references{
Morrison, Julie L., et al. "GeneRank: using search engine technology for the analysis of microarray experiments." BMC bioinformatics 6.1 (2005): 233.\cr
Page, Lawrence, et al. "The PageRank citation ranking: bringing order to the web." (1999).\cr
}
\author{
Yupeng Cun \email{yupeng.cun@gmail.com}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See Also as \code{classify.frsvm}
}
\examples{
#See Also as {classify.frsvm}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ FrSVM  }
%%\keyword{ GeneRank }% __ONLY ONE__ keyword per line
