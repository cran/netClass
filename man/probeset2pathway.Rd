\name{probeset2pathway}
\alias{probeset2pathway}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generae a mean gene expression of genes of each pathway matrix
}
\description{
Generae a mean gene expression of genes of each pathway matrix
}
\usage{
probeset2pathway(x = x, int = int, sigGens = sigGens)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
gene expression data
}
  \item{int}{
common genes between pathway genes and genes in gene expression profile
}
  \item{sigGens}{
significant gene expression using SAM methods
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{kse }{an matrix with n pathways and p samples}
%%  \item{comp2 }{Description of 'comp2'}
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
See Also as \code{classify.aep}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ aep }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
