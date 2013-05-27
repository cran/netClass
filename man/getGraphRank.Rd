\name{getGraphRank}
\alias{getGraphRank}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Random walk kernel matrix smoothing t-statistic
}
\description{
Using Random walk kernel matrix of network to smooth t-statistic of each gene
}
\usage{
getGraphRank(x = x, y = y, Gsub = Gsub, sca = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a matrix of expression measurements with p samples and n genes.
}
  \item{y}{
a factor of length p comprising the class labels. 
}
  \item{Gsub}{
 Random Walk Kernel matrix of network
}
  \item{sca}{
Sacling data or not
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{r }{return a smoothed t-statistic of each gene'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Yupeng Cun, Holger Frohlich (2013) Network and Data Integration for Biomarker Signature Disvovery via Network Smoothed T-Statistics
}
\author{
Yupeng Cun \email{yupeng.cun@gmail.com}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See Also as \code{getGraphRank}
}
\examples{
#See also \code{classfy.stsvm}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ stSVM }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
