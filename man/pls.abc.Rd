\name{pls.abc}
\alias{pls.abc}
\alias{Routput_test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ABC inference using PLS components computed from data.
}
\description{
The function performs ABC inference using PLS transformed summary statistics.
}
\usage{
pls.abc(obs, param, sumstats, obspar=NULL, abcmethod=abc,transfile = "Routput_test", 
bc=FALSE, err.only=TRUE, errfn=rsse,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obs}{(matrix of) observed summary statistics.}
  \item{param}{matrix of simulated model parameter values.}
  \item{sumstats}{matrix of simulated summary statistics.}
  \item{obspar}{optional observed parameters (for use to assess simulation performance).}
  \item{abcmethod}{a function to perform ABC inference, e.g. the \code{abc} function from package \emph{abc}.}
  \item{transfile}{path to file containing the PLS transformation (see Details section and documentation for \emph{ABCtoolbox}).}
  \item{bc}{a boolean value indicating whether the Box-Cox transformation should be applied to the statistics prior to transformation.}
  \item{err.only}{a boolean value indicating whether only the simulation error should be returned.  Note: if \code{err.only=TRUE}, \code{obspar} must be supplied.}
  \item{errfn}{an error function to assess ABC inference performance.}
  \item{\dots}{any other optional arguments to the ABC inference procedure (e.g. arguments to the \code{abc} function).}
}
\details{The function uses a precomputed PLS transformation file to construct new summary statistics, being linear combinations of the original summary statistics.  This is achieved using the \code{transformer} command line script supplied with the \emph{ABCtoolbox} software package.  Note that the precomputed transformation file can be created for example, with the \emph{find_pls.r} R script supplied in \emph{ABCtoolbox}.  An example \code{Routput_test} for three components is supplied as a dataset.  See the package documentation for more information on these two script files.}
\value{
A list with the following components:
  \item{post.sample}{an array of dimension \code{nacc x npar x ndatasets} giving the posterior sample for each observed dataset.  Not returned if \code{err.only=TRUE}.}
  \item{err}{simulation error (if \code{obspar} is supplied).}
}
\references{
Blum, M. G. B, Nunes, M. A., Prangle, D. and Sisson, S. A. (2013) A comparative review of dimension reduction methods in approximate Bayesian 
computation. \emph{Stat. Sci. (to appear)}.\cr\cr
Wegmann, D. et al. (2010) ABCtoolbox: A versatile toolkit for approximate Bayesian computations. \emph{BMC Bioinformatics} \bold{11}, 116--123.\cr\cr
Wegmann, D. et al. (2009) Efficient approximate Bayesian computation coupled With Markov chain Monte Carlo Without Likelihood. \emph{Genetics} \bold{182} (4), 1207--1218.
}
\author{Matt Nunes}
\note{This function requires that the \code{transformer} script from ABCtoolbox is in the system path.}

%\section{Warning}{} 

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{

# example transformation file:

data(Routput_test)
write.table(Routput_test,file="Routput_test",quote=FALSE,row.names=FALSE,col.names=FALSE)

data(coal)

\dontrun{
pls.abc(coal[1:3,3:7], coal[,1:2], coal[,3:7], transfile = "Routput_test",tol=.1,
method="rejection")
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
