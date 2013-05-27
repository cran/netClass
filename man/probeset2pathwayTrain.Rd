\name{probeset2pathwayTrain}
\alias{probeset2pathwayTrain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Search CROG in training data
}
\description{
Search CROG in training data, and using these CORG set to make a matrix for pathways.

}
\usage{
probeset2pathwayTrain(x = x, y = y, int = int)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
gene expression data
}
  \item{y}{
a factor of length p comprising the class labels.
}
  \item{int}{
Common genes between gene expression data and interaction network.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{ap }{top ranked pathays}
\item{selectedGenes }{ CROG genes }
....
}
\references{
Lee E, Chuang H-Y, Kim J-W, Ideker T, Lee D (2008) Inferring Pathway Activity toward Precise Disease Classification. PLoS Comput Biol 4(11): e1000217. doi:10.1371/journal.pcbi.1000217
}
\author{
Yupeng Cun \email{yupeng.cun@gmail.com}

}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See Also as \code{pac.cv}
}
\examples{
#See Also as \name{pac.cv}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ pac }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
