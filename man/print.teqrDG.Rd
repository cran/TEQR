\name{print.teqrDG}
\alias{print.teqrDG}
\title{print.teqrDG}
\description{Print file for the class teqrDG.}
\usage{
\method{print}{teqrDG}(x,\dots)
}
\arguments{
\item{x}{object of class teqrDG}
\item{...}{additional arguments to be passed to lower level functions}
}
\value{The output from print.teqrDG is the dose escalation/expansion/de-escalation guidelines table.}
\author{M.S. Blanchard <sblanchard@coh.org>}
\examples{
output<-teqrDG(TotalN=12,pTarget=0.2,eq1=0.05,eq2=0.05,tootoxic=0.34)
output
}
\keyword{print} 


