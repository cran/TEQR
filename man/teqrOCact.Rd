\name{teqrOCact}
\alias{teqrOCact}
\title{TEQR activity design trial simulator}
\description{
     teqrOCact: Runs clinical trial simulations for the TEQR Design implemented as a phase I activity design constrained for toxicity and provides a listing of the operating characteristics.
     
     The target equivalence range (TEQR) design, previously named toxicity equivalence range design is an up-down design patterned after the modified toxicity probability interval design. The TEQR design extends the well known 3+3 design to allow for: an explicit target range for the dose limiting toxicity (DLT) rate,  more than 6 subjects at the maximum tolerated dose (MTD), and specification of a too-toxic rate, which closes a dose level. In addition, the TEQR design can be implemented to target a dose based either on toxicity rate or activity rate constrained for toxicity.
     
     The TEQR is better than the 3+3 when compared on: 1) the number of times the dose at or nearest the target toxicity level was selected as the MTD, 2) the number of subjects assigned to doses levels, at, or nearest the MTD, and 3) the overall trial DLT rate. For sample sizes at the MTD greater that 6 TEQR more accurately and more precisely estimates the rate of toxicity at the MTD because a larger number of subjects are studied at the MTD dose. The TEQR on average uses fewer subjects and provide reasonably comparable results to the  continual reassessment method (CRM) in the number of times the dose at or nearest the target toxicity level was selected as the MTD and the number of subjects assigned doses, at, or nearest that target and in overall DLT rate. 
}
\usage{
teqrOCact(sim,firstdose,proba,probc,cohortSize,MaxNoCohorts,
RP2Dss,pTarget,eq1,eq2,tootoxic,toxcon,llactivity)}
\arguments{
\item{sim}{Number of trial simulations desired.}
\item{firstdose}{Dose level where the trial will start. The default setting is 2.}
\item{proba}{Vector of activities associated with the dose elvels from 1 to n}
\item{probc}{Vector of toxicities associated with the dose levels from 1 to n 
(e.g probt=c(.05,.15,.23,.34,.51.76) for dose levels 1 to 6.)}
\item{cohortSize}{Number of subjects per cohort. The default setting is 3.}
\item{MaxNoCohorts}{Highest number of cohorts possible. The default setting is 30.}
\item{RP2Dss}{Desired sample size at the RP2D. The default setting is 12.} 
\item{pTarget}{Target activity probability at the MTD.}
\item{eq1}{Range subtracted from pTarget to achieve the lower boundary for the equivalence range.}
\item{eq2}{Range added to pTarget to achieve the upper boundary for the equivalence range.}
\item{tootoxic}{Toxicity probability deemed too high to allow further study at a dose level. Default is set to 1.01}
\item{toxcon}{toxicity constaint level}
\item{llactivity}{Required lower limit of activity to end the trial}
}
\value{teqrOCact outputs the trial simulation operating characteristics, the simulation data sets, and all the simulation parameters.  The simulation operating characteristics include, i) the average number of subjects studied at each dose level, ii) the rate the dose level is chosen as the MTD (including if there is no MTD identified), iii) median study sample size, iv) mean study DLT rate, v) mean DLT rate at the MTD, vi) mean 95\% confidence interval rate at the MTD, vii) proportion of the trials with sample size at or above the level desired at the MTD, viii) number of simulated trials that do not determine an MTD and ix) number of simulated trials.
}
\references{Blanchard MS and Longmate JA, (2010). Toxicity equivalence range design (TEQR): A practical Phase I design. \emph{Contemp Clin Trials} (2010), doi:10.1016/j.cct.2010.09.011}
\author{M.S. Blanchard <sblanchard@coh.org>}
\examples{
set.seed(1111)
OperChar<-teqrOCact(
sim=500,
firstdose=2,
proba=c(.10,.65,.85), 
probc=c(.10,.20,.30),
cohortSize=3,
MaxNoCohorts=30,
RP2Dss=12,
pTarget=.75,  
eq1=.10,
eq2=.25,
toxcon=.34,
llactivity=.30)
OperChar
}
\keyword{design} 
