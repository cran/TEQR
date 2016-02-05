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
probt=c(0.05,0.15,0.23,0.34,0.51,0.76), 
cohortSize=3,
MaxNoCohorts=30,
MTDss=12,
pTarget=0.2,
eq1=0.05,
eq2=0.05,
tootoxic=0.34)
OperChar
}
\keyword{print} 


