\name{teqr.MTDdatatox}
\alias{teqr.MTDdatatox}
\title{teqrOCtox intermediate file that creates MTD level data}
\description{teqr.MTDdatatox takes the simData and DLdata output and creates the MTD level data.
}
\usage{
teqr.MTDdatatox(simData = simData, DLdata = DLdata)
}
\arguments{
  \item{simData}{includes the simulation results data,simresult, cohortSize,probt, MTDss, pTarget, lowerlimit, upperlimit,tootoxic, the number of simulation (sim)}
  \item{DLdata}{simNo,doselevel,stox (no of subjects experiencing a toxicity),dllength (length of the cohort),toxl (lower 95\% CI),toxu (upper 95\% CI) ,toxest (toxicity estimate)}
}
\references{Blanchard MS and Longmate JA, (2011). Toxicity equivalence range design (TEQR): A practical Phase I design. \emph{Contemporary Clinical Trials,(32)}, 114-121.}
\author{M.S. Blanchard <sblanchard@coh.org>}
\keyword{design}