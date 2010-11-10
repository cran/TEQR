\name{teqrOC}
\alias{teqrOC}
\title{TEQR Design trial simulator}
\description{
     teqrOC: Runs clinical trial simulations for the TEQR Design and provides a listing of the operating characteristics.

     The target equivalence range (TEQR) design is a frequentist implementation of the modified toxicity probability interval (mTPI) design and a competitor to the standard 3+3 design (3+3).  The 3+3 is the work horse design in Phase I. It is good at determining if a safe dose exits, but provides poor accuracy and precision in estimating the level of toxicity at the maximum tolerated dose (MTD).  The TEQR is better than the 3+3 when compared on: 1) the number of times the dose at or nearest the target toxicity level was selected as the MTD, 2) the number of subjects assigned to doses levels, at, or nearest the MTD, and 3) the overall trial DLT rate. TEQR more accurately and more precisely estimates the rate of toxicity at the MTD because a larger number of subjects are studied at the MTD dose. The TEQR on average uses fewer subjects and provide reasonably comparable results to the  continual reassessment method (CRM) in the number of times the dose at or nearest the target toxicity level was selected as the MTD and the number of subjects assigned doses, at, or nearest that target and in overall DLT rate. 
}
\usage{
     teqrOC(sim,firstdose,probt,cohortSize,MaxNoCohorts,MTDss,pTarget,eq1,eq2,tootoxic)}
\arguments{
\item{sim}{Number of trial simulations desired.}
\item{firstdose}{Dose level where the trial will start. The default setting is 2.}
\item{probt}{Vector of toxicities associated with the dose levels from 1 to n 
(e.g probt=c(.05,.15,.23,.34,.51.76) for dose levels 1 to 6.). We recommend that there be at least 4 dose levels with at least one level below and one level above the level expected to be at the target toxicitiy probability.}
\item{cohortSize}{Number of subjects per cohort. The default setting is 3.}
\item{MaxNoCohorts}{Highest number of cohorts possible. The default setting is 30.}
\item{MTDss}{Desired sample size at the MTD. The default setting is 12.} 
\item{pTarget}{Target toxicity probability at the MTD.}
\item{eq1}{Range subtracted from pTarget to achieve the lower boundary for the equivalence range.}
\item{eq2}{Range added to pTarget to achieve the upper boundary for the equivalence range.}
\item{tootoxic}{Toxicity probability deemed too high to allow further study at a dose level.}
}
\value{teqrOC outputs the trial simulation operating characteristics, the simulation data sets, and all the simulation parameters.  The simulation operating characteristics include, i) the average number of subjects studied at each dose level, ii) the rate the dose level is chosen as the MTD (including if there is no MTD identified), iii) median study sample size, iv) mean study DLT rate, v) mean DLT rate at the MTD, vi) mean 95\% confidence interval rate at the MTD, vii) proportion of the trials with sample size at or above the level desired at the MTD, viii) number of simulated trials that do not determine an MTD and ix) number of simulated trials.
}
\references{Blanchard MS and Longmate JA, (2010). Toxicity equivalence range design (TEQR): A practical Phase I design. \emph{Contemp Clin Trials} (2010), doi:10.1016/j.cct.2010.09.011}
\author{M.S. Blanchard <sblanchard@coh.org>}
\examples{
set.seed(1111)OperChar<-teqrOC(
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
\keyword{design} 

