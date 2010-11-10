\name{print.teqrOC}
\alias{print.teqrOC}
\title{print.teqrOC}
\description{Print file for the class teqrOC.}
\usage{
\method{print}{teqrOC}(x,\dots)
}
\arguments{
\item{x}{object of class teqrOC}
\item{...}{additional arguments to be passed to lower level functions}
}
\value{The output from print.teqr is a listing of the operating characteristics for the design as specified by the user.}
\author{M.S. Blanchard <sblanchard@coh.org>}
\examples{
set.seed(1111)
OperChar<-teqrOC(
sim=1000,
firstdose=2,
probt=c(.05,.15,.23,.34,.51,.76), 
cohortSize=3,
MaxNoCohorts=30,
MTDss=18,
pTarget=.2,
eq1=.05,
eq2=.05,
tootoxic=.34)
OperChar
}
\keyword{print} 


