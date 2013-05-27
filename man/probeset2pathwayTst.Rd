\name{probeset2pathwayTst}
\alias{probeset2pathwayTst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Applied CROG to testing data
}
\description{
Applied CORG and pathways activities lists to make a matrix for pathways for test data.
}
\usage{
probeset2pathwayTst(x = x, apTrain = apTrain)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
gene expression data
}
  \item{apTrain}{
PAC objext which contain CORG and pathways activities lists of training data.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{ap }{top ranked pathays}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
See Also as \code{pac.cv, probeset2pathwayTrain}
}
\examples{
#See Also as \code{pac.cv, probeset2pathwayTrain}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ pac }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
