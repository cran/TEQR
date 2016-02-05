\name{teqr.MTDdataconACT}
\alias{teqr.MTDdataconACT}
\title{teqrOCact intermediate file that creates MTD/RP2D level data}
\description{teqr.MTDdatatox takes the simData and DLdata output and creates the MTD level data.}
\usage{
teqr.MTDdataconACT(simData = simData, DLdata = DLdata)
}
\arguments{
  \item{simData}{includes the simulation results data, cohortSize,proba,probc, MTDss,pTarget,lowerlimit,upperlimit,tootoxic, toxcon, llactivity, number of simulations (sim)}
  \item{DLdata}{simNo,doselevel,stox (no of subjects experiencing activity),dllength (length of the cohort),toxl (lower 95\% CI),toxu (upper 95\% CI) ,toxCest (toxicity estimate)}
}
\references{Blanchard MS and Longmate JA, (2011). Toxicity equivalence range design (TEQR): A practical Phase I design. \emph{Contemporary Clinical Trials,(32)}, 114-121.}
\author{M.S. Blanchard <sblanchard@coh.org>}
\keyword{design}

