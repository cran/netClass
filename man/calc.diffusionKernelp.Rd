\name{calc.diffusionKernelp}
\alias{calc.diffusionKernelp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computing the Random Walk Kernel matrix of network
}
\description{
 Computing the Random Walk Kernel matrix of network 
}
\usage{
calc.diffusionKernelp(L, is.adjacency = TRUE, p = 3, a = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{L}{
an adjacency matrix that represents the underlying biological network.
}
  \item{is.adjacency}{
using adjacency of graph or not
}
  \item{p}{
#(p) random walk step(s) of random walk kernel
}
  \item{a}{
constant value of random walk kernel
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  \item{R}{Return a Random Walk Kernel matrix of given network, L.}

}
\references{
Kondor, R. I., & Lafferty, J. (2002, July). Diffusion kernels on graphs and other discrete input spaces. In MACHINE LEARNING-INTERNATIONAL WORKSHOP THEN CONFERENCE- (pp. 315-322).
}
\author{
Yupeng Cun \email{yupeng.cun@gmail.com}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See Also as \code{classify.stsvm}
}
\examples{
library(netClass)
data(ad.matrix) 
#dk= calc.diffusionKernelp(L=ad.matrix, is.adjacency=TRUE, p=2,a=1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{stSVM }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
