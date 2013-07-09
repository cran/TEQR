\name{desc}
\alias{desc}
\title{Descriptive Statistics
}
\description{
  desc caculates select descriptive statistics for the variable y.
}
\usage{
desc(y, pcts=c(0.025,0.05,0.95,0.975), nsig=4)
}
\arguments{
\item{y}{variable of interest}
\item{pcts}{selected percentiles, the default is c(0.025,0.05,0.95,0.975) }
\item{nsig}{number of significant figures}
}

\value{desc prints descriptive statistics including mean, median, standard deviation,range, min, max and selected percentiles.
}

\author{
M.S. Blanchard <sblanchard@coh.org>
}

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
TotalN<-rep(NA,500)
for (i in 1:500){
TotalN[i]<-sum(OperChar$DLdata$dllength[OperChar$DLdata$simNo==i])
}
desc(TotalN)
}


\keyword{statistics}
