\name{train.stsvm}
\alias{train.stsvm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Training the data using stsvm methods
}
\description{
Training the data using stsvm methods
}
\usage{
train.stsvm(x=x, y=y, DEBUG=FALSE,Gsub=Gsub, op.method="sp", op=10,aa=100,
			dk=dk, dk.tf=0.05,seed = 1234,Cs=10^(-3:3),EN2SY=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
expression data for training
}
  \item{y}{
Class labels
}
  \item{DEBUG}{
show debugging information in screen more or less.
}
  \item{Gsub}{
 an adjacency matrix that represents the underlying biological network.
}
\item{op.method}{
Method for selecet optimal feature subgoups: pt is permutation test, sp is span bound.
}
  \item{op}{
optimal on top op% ranked genes
}
  \item{aa}{
permutation test steps
}
  \item{dk}{
Random Walk Kernel matrix of network
}
  \item{dk.tf}{
cut off p-value of permutation test
}
  \item{seed}{
seed for random sampling.
}
  \item{Cs}{
Soft-margin tuning parameter of the SVM. Defaults to \code{10^c(-3:3)}.
}
 \item{EN2SY}{
A list for mapping gene sybol ids or entez ids.
}

}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
The list returned
%%  If it is a LIST, use
\item{trained }{The tranined models for traning folds}
\item{feat }{The feature selected by each by the train }
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
See cv.stsvm
}
\examples{
#see cv.stsvm

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ stsvm }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line%
