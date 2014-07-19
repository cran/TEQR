\name{teqrOCtox}
\alias{teqrOCtox}
\title{TEQR Design trial simulator}
\description{
     teqrOCtox: Runs clinical trial simulations for the TEQR Design implemented as a Phase I toxicity design and provides a listing of the operating characteristics.

The target equivalence range (TEQR) design, previously named toxicity equivalence range design is an up-down design patterned after the modified toxicity probability interval design. The TEQR design extends the well known 3+3 design to allow for: an explicit target range for the dose limiting toxicity (DLT) rate,  more than 6 subjects at the maximum tolerated dose (MTD), and specification of a too-toxic rate, which closes a dose level. In addition, the TEQR design can be implemented to target a dose based either on toxicity rate or activity rate constrained for toxicity.

The TEQR is better than the 3+3 when compared on: 1) the number of times the dose at or nearest the target toxicity level was selected as the MTD, 2) the number of subjects assigned to doses levels, at, or nearest the MTD, and 3) the overall trial DLT rate. For sample sizes at the MTD greater that 6 TEQR more accurately and more precisely estimates the rate of toxicity at the MTD because a larger number of subjects are studied at the MTD dose. The TEQR on average uses fewer subjects and provide reasonably comparable results to the  continual reassessment method (CRM) in the number of times the dose at or nearest the target toxicity level was selected as the MTD and the number of subjects assigned doses, at, or nearest that target and in overall DLT rate. 


}
\usage{
     teqrOCtox(sim,firstdose,probt,cohortSize,MaxNoCohorts,MTDss,pTarget,eq1,eq2,tootoxic)}
\arguments{
\item{sim}{Number of trial simulations desired. We recommend at least 1000. The example runs 500 to meet R package time limits.}
\item{firstdose}{Dose level where the trial will start. The default setting is 2.}
\item{probt}{Vector of toxicities associated with the dose levels from 1 to n 
(e.g probt=c(.05,.15,.23,.34,.51,.76) for dose levels 1 to 6.). We recommend that there be at least 4 dose levels with at least one level below and one level above the level expected to be at the target toxicitiy probability.}
\item{cohortSize}{Number of subjects per cohort. The default setting is 3.}
\item{MaxNoCohorts}{Highest number of cohorts possible. The default setting is 30.}
\item{MTDss}{Desired sample size at the MTD. The default setting is 12.} 
\item{pTarget}{Target toxicity probability at the MTD.}
\item{eq1}{Range subtracted from pTarget to achieve the lower boundary for the equivalence range.}
\item{eq2}{Range added to pTarget to achieve the upper boundary for the equivalence range.}
\item{tootoxic}{Toxicity probability deemed too high to allow further study at a dose level.}
}
\value{teqrOC outputs the trial simulation operating characteristics, the simulation data sets, and all the simulation parameters.  The simulation operating characteristics include, i) the average number of subjects studied at each dose level, ii) the rate the dose level is chosen as the MTD (including if there is no MTD identified), iii) median study sample size, iv) mean study DLT rate, v) mean DLT rate at the MTD, vi) mean 95\% confidence interval rate at the MTD, vii) proportion of the trials with sample size at or above the level desired at the MTD, viii) number of simulated trials that do not determine an MTD and ix) number of simulated trials.

Datasets:

simData$simresults gives a line of data for each cohort.
The columns are:
simNo - simulation number, 
doselevel - dose level, 
tox - the number of subjects the experience a toxicity in the cohort, and 
cumtox - cumulative toxicity in the cohort.

DLdata gives a line of data for each dose level. The columns are:
simNo - simulation number,
doselevel - dose level,
stox - number of subjects that experience a toxicity at that doselevel,
dllength - number of subjects studied at a doselevel,
toxl - lower limit of the 95\% confidence limits for toxicity rate,
toxu - upperlimit of the 95\% confidence limits for toxicity rate, and
toxest - toxicity rate estimate at that doselevel.

MTDdata gives a line of data at the MTD dose selected. The columns names are the same as for DLdata.
}
\references{Blanchard MS and Longmate JA, (2011). Toxicity equivalence range design (TEQR): A practical Phase I design. \emph{Contemporary Clinical Trials,(32)}, 114-121.}
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
\keyword{design} 

