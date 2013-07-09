\name{teqrDG}
\alias{teqrDG}
\title{TEQR Design Dose Escalation/Expansion/De-escalation Guidelines}
\description{
     teqrDG prints the TEQR Design Dose Escalation/Expansion/De-escalation Guidelines.

The target equivalence range (TEQR) design, previously named toxicity equivalence range design is an up-down design patterned after the modified toxicity probability interval design. The TEQR design extends the well known 3+3 design to allow for: an explicit target range for the dose limiting toxicity (DLT) rate,  more than 6 subjects at the maximum tolerated dose (MTD), and specification of a too-toxic rate, which closes a dose level. In addition, the TEQR design can be implemented to target a dose based either on toxicity rate or activity rate constrained for toxicity.

The TEQR is better than the 3+3 when compared on: 1) the number of times the dose at or nearest the target toxicity level was selected as the MTD, 2) the number of subjects assigned to doses levels, at, or nearest the MTD, and 3) the overall trial DLT rate. For sample sizes at the MTD greater that 6 TEQR more accurately and more precisely estimates the rate of toxicity at the MTD because a larger number of subjects are studied at the MTD dose. The TEQR on average uses fewer subjects and provide reasonably comparable results to the  continual reassessment method (CRM) in the number of times the dose at or nearest the target toxicity level was selected as the MTD and the number of subjects assigned doses, at, or nearest that target and in overall DLT rate. 

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

