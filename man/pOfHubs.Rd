\name{pOfHubs}
\alias{pOfHubs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computing p value of hubs using the permutation test 
}
\description{
Computing p value of hubs using the permutation test 
}
\usage{
pOfHubs(x = x, y = y, gHub = gHub, hubs = hubs, nperm = nperm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
gene expression data
}
  \item{y}{
a factor of length p comprising the class labels.
}
  \item{gHub}{
Subgraph of hubs of graph Gs
}
  \item{hubs}{
Hubs in graph Gs
}
  \item{nperm}{
number of permutation test steps
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{pVal}{Permutation test Pvalues of each hub}
  \item{hub }{name of hubs}
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

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
# see \code{pOfHubs}


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hubc }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
