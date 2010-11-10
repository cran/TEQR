\name{teqrDG}
\alias{teqrDG}
\title{TEQR Design Dose Escalation/Expansion/De-escalation Guidelines}
\description{
     teqrDG prints the TEQR Design Dose Escalation/Expansion/De-escalation Guidelines.

     The target equivalence range (TEQR) design is a frequentist implementation of the modified toxicity probability interval (mTPI) design and a competitor to the standard 3+3 design (3+3).  The 3+3 is the work horse design in Phase I. It is good at determining if a safe dose exits, but provides poor accuracy and precision in estimating the level of toxicity at the maximum tolerated dose (MTD).  The TEQR is better than the 3+3 when compared on: 1) the number of times the dose at or nearest the target toxicity level was selected as the MTD, 2) the number of subjects assigned to doses levels, at or nearest the MTD, and 3) the overall trial DLT rate. TEQR more accurately and more precisely estimates the rate of toxicity at the MTD because a larger number of subjects are studied at the MTD dose. The TEQR on average uses fewer subjects and provide reasonably comparable results to the  continual reassessment method (CRM) in the number of times the dose at or nearest the target toxicity level was selected as the MTD and the number of subjects assigned doses, at, or nearest the target and in overall DLT rate. 
}
\usage{
     teqrDG(TotalN,pTarget,eq1,eq2,tootoxic)}
\arguments{
\item{TotalN}{Total Number of subjects being considered at a dose level}
\item{pTarget}{Target toxicity probability at the MTD.}
\item{eq1}{Range subtracted from pTarget to achieve the lower boundary for the equivalence range.}
\item{eq2}{Range added to pTarget to achieve the upper boundary for the equivalence range.}
\item{tootoxic}{Probability cutoff at which the toxicity probability is deemed too high to allow further study at that dose level.}
}
\value{TEQR outputs a table of dose escalation guidelines (DoseGuideTable) and a table of the underlying probabilities (ProbTable).}
\references{Blanchard MS and Longmate JA, (2010). Toxicity equivalence range design (TEQR): A practical Phase I design. \emph{Contemp Clin Trials} (2010), doi:10.1016/j.cct.2010.09.011}
\author{M.S. Blanchard <sblanchard@coh.org>}
\examples{
output<-teqrDG(TotalN=12,pTarget=.2,eq1=.05,eq2=.05,tootoxic=.34)
output
}
\keyword{design} 

