\name{print.teqrOCact}
\alias{print.teqrOCact}
\title{print.teqrOCact}
\description{Print file for the class teqrOCact.}
\usage{
\method{print}{teqrOCact}(x,\dots)
}
\arguments{
\item{x}{object of class teqrOCact}
\item{...}{additional arguments to be passed to lower level functions}
}
\value{The output from print.teqr is a listing of the operating characteristics for the design as specified by the user.}
\author{M.S. Blanchard <sblanchard@coh.org>}
\examples{
set.seed(1111)
OperChar<-teqrOCact(
sim=1000,
firstdose=2,
proba=c(0.50,0.80,0.90), 
probc=c(0.10,0.10,0.10),
cohortSize=3,
MaxNoCohorts=30,
RP2Dss=6,
pTarget=1.00,  
eq1=0.34,
eq2=0.0,
toxcon=0.34,
llactivity=0.66)
OperChar
}
\keyword{print} 
