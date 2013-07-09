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
proba=c(.50,.80,.90), 
probc=c(.10,.10,.10),
cohortSize=3,
MaxNoCohorts=30,
RP2Dss=6,
pTarget=.75,  
eq1=.10,
eq2=.25,
toxcon=.34,
llactivity=.30)
OperChar
}
\keyword{print} 
