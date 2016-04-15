\name{deconvolveComp}
\alias{deconvolveComp}
\title{Deconvolution of compounds in samples}
\description{Deconvolution of GC-MS data}
\usage{
deconvolveComp(Experiment, decParameters, 
samples.to.process=NULL)
}
\arguments{
  \item{Experiment}{A 'MetaboSet' S4 object containing the experiment data previously created by newExp.}
  \item{decParameters}{The software deconvolution parameters object previously created by setDecPar}
  \item{samples.to.process}{Vector indicating which samples are to be processed.} 
}
\details{
See eRah vignette for more details. To open the vignette, execute the following code in R:
vignette("eRahManual", package="erah")

(Available also in http://metabolomicsplatform.com/applications).
}
\value{
The function returns an updated S4 'MetaboSet' class, where the GC-MS samples have been now deconvolved.
}
\references{
[1] eRah: an R package for spectral deconvolution, alignment, and metabolite identification in GC/MS-based untargeted metabolomics. Xavier Domingo-Almenara, Alexandre Perera, Maria Vinaixa, Sara Samino, Xavier Correig, Jesus Brezmes, Oscar Yanes. (2016) Article in Press.
}
\author{
Xavier Domingo-Almenara. xavier.domingo@urv.cat
}

\seealso{
\code{\link{newExp}}, \code{\link{setAlPar}}
}
\examples{

# Deconvolve data from a created experiment by \code{\link{newExp}}.
# ex <- newExp(instrumental="path")

# The following will set eRah for analyzing the chromatograms
# from minutes 5 to 15, and withouth taking into account the masses
# 35:69,73:75,147:149, widht a minimum peak widht of 0.7 seconds.

ex.dec.par <- setDecPar(min.peak.width=0.7, min.peak.height=5000, 
noise.threshold=500, avoid.processing.mz=c(35:69,73:75,147:149), 
analysis.time=c(5,15))

# An now deconvolve the compounds in the samples:
# ex <- deconvolveComp(ex, decParameters=ex.dec.par)

}


