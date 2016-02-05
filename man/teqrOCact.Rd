\name{teqrOCact}
\alias{teqrOCact}
\title{Activity Constrained for Toxicity (ACT) design trial simulator}
\description{
The activity constrained for toxicity (ACT) design implements a cumulative cohort design in conjunction with safety rules similar to those from the toxicity equivalence range design or the modified toxicity probability interval design. The unique feature of this design is that dose is escalated based on lack of activity rather than on lack of toxicity and is de-escalated only if an unacceptable level of toxicity is experienced.     
}
\usage{
teqrOCact(sim,firstdose,proba,probc,cohortSize,MaxNoCohorts,
RP2Dss,pTarget,eq1,eq2,tootoxic,toxcon,llactivity)}
\arguments{
\item{sim}{Number of trial simulations desired. We recommend at least 1000. The example runs 500 to meet R package time limits.}
\item{firstdose}{Dose level where the trial will start. The default setting is 2.}
\item{proba}{Vector of activities associated with the dose elvels from 1 to n}
\item{probc}{Vector of toxicities associated with the dose levels from 1 to n 
(e.g probt=c(0.05,0.15,0.23,0.34,0.51,0.76) for dose levels 1 to 6.)}
\item{cohortSize}{Number of subjects per cohort. The default setting is 3.}
\item{MaxNoCohorts}{Highest number of cohorts possible. The default setting is 30.}
\item{RP2Dss}{Desired sample size at the recommended phase 2 dose (RP2D). The default setting is 12.} 
\item{pTarget}{Target activity probability at the RP2D.}
\item{eq1}{Range subtracted from pTarget to achieve the lower boundary for the equivalence range.}
\item{eq2}{Range added to pTarget to achieve the upper boundary for the equivalence range.}
\item{tootoxic}{Toxicity probability deemed too high to allow further study at a dose level. Default is set to 1.01 and should not be changed in studing the ACT design.}
\item{toxcon}{toxicity constaint level}
\item{llactivity}{Required lower limit of activity to end the trial}
}
\value{teqrOCact outputs the trial simulation operating characteristics, the simulation data sets, and all the simulation parameters.  The simulation operating characteristics include, i) the average number of subjects studied at each dose level, ii) the rate the dose level is chosen as the RP2D (including if there is no RP2D identified), iii) median study sample size, iv) mean study event rate, v) mean event rate at the RP2D, vi) mean 95\% confidence interval length for the event rate at the RP2D, vii) mean toxicity rate at the the RP2D, viii) proportion of the trials with sample size at or above the level desired at the RP2D, viii) number of simulated trials that do not determine an RP2D and ix) number of simulated trials.

Datasets:

simData$simresults gives a line of data for each cohort.
The columns are:
simNo - simulation number, 
doselevel - dose level, 
act - the number of subjects that achieve the activity endpoint in the cohort,
toxC - the number of subjects that experience a toxicity in the cohort,
cumact - the cumulative rate of subjects achieving the activity endpoint at that dose level, and
cumtoxC - the cumulative rate of subjects experiencing toxicity at that dose level.

DLdata gives a line of data for each dose level. The columns are:
simNo - simulation number,
doselevel - dose level,
sact - sum of the number of subjects that achieve the activity endpoint at that doselevel,
stoxC - sum of the number of subjects that experience the toxicity endpoint at that dose level,  
dllength - number of subjects studied at a doselevel,
actl - lower limit of the 95\% confidence limits for the activity rate,
actu - upperlimit of the 95\% confidence limits for the activity rate, 
actest - activity rate estimate at that doselevel, and
toxCests - toxicity rate estimate at that doselevel.

RP2Ddata gives a line of data at the RP2D dose selected. The columns names are the same as for DLdata.
}

\author{M.S. Blanchard <sblanchard@coh.org>}
\examples{
set.seed(1111)
OperChar<-teqrOCact(
sim=500,
firstdose=2,
proba=c(0.65,0.75,0.85,1.00), 
probc=c(0.10,0.20,0.30,0.40),
cohortSize=3,
MaxNoCohorts=30,
RP2Dss=12,
pTarget=1.00,  
eq1=0.34,
eq2=0.0,
toxcon=0.34,
llactivity=0.66)
OperChar
}
\keyword{design} 
