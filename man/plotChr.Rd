\name{plotChr}
\alias{plotChr}
\title{Plotting sample chromatogram}
\description{Plot the sample chromatogram.}
\usage{
plotChr(Experiment, N.sample=1, type=c("BIC","TIC","EIC"), 
xlim=NULL, mz=NULL)
}
\arguments{
  \item{Experiment}{A 'MetaboSet' S4 object containing the experiment.}
  \item{N.sample}{Integer. The number of the sample to query.}
  \item{type}{The type of plotting, Base Ion Chromatogram (BIC), Total Ion Chromatogram (TIC), or Extracted Ion Chromatogram (EIC).}
  \item{xlim}{The range in minutes, separated by comas: c(rt.min, rt.max) of the limits of plotting. By default, all the chromatogram is plotted.}
  \item{mz}{Just when EIC is selected. The range separated by comas: c(mz.min, mz.max) or a vector of numbers: c(50,67,80), of the masses to be ploted.}
}

\examples{
\dontrun{
# First, an experiment has to be already created by newExp()
# then, each sample chromatogram can be plotted by:

plotChr(Experiment, 1, "BIC")
plotChr(Experiment, 1, "TIC", xlim=c(5,7))  #Plots from minute 5 to 7.

plotChr(Experiment, 1, "EIC", mz=50:70, xlim=c(5,7))  #Plots 
# from minute 5 to 7, and only the masses from 50 to 70.

plotChr(Experiment, 1, "EIC", xlim=c(7,7.5), mz=c(50,54,70))  #Plots
# the EIC from minute 7 to 7.5, and only the masses 50, 54 and 70.

}}

\seealso{
\code{\link{sampleInfo}}

}