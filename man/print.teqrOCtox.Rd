\name{print.teqrOCtox}
\alias{print.teqrOCtox}
\title{print.teqrOCtox}
\description{Print file for the class teqrOCtox.}
\usage{
\method{print}{teqrOCtox}(x,\dots)
}
\arguments{
\item{x}{object of class teqrOCtox}
\item{...}{additional arguments to be passed to lower level functions}
}
\value{The output from print.teqrOCtox is a listing of the operating characteristics for the design as specified by the user.}
\author{M.S. Blanchard <sblanchard@coh.org>}
\examples{
set.seed(1111)

OperChar<-teqrOCtox(
sim=500,
firstdose=2,
probt=c(.05,.15,.23,.34,.51,.76), 
cohortSize=3,
MaxNoCohorts=30,
MTDss=12,
pTarget=.2,
eq1=.05,
eq2=.05,
tootoxic=.34)
OperChar
}
\keyword{print} 


