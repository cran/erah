\name{alignComp}
\alias{alignComp}
\title{Alignment of compounds}
\description{Alignment of GC-MS deconvolved compounds}
\usage{
alignComp(Experiment, alParameters, blocks.size=NULL)
}
\arguments{
  \item{Experiment}{A 'MetaboSet' S4 object containing the experiment data previously created by newExp and deconvolved by deconvolveComp.}
  \item{alParameters}{The software alignment parameters object previously created by setAlPar}
  \item{blocks.size}{For experiment of more than 1000 samples, and depending on the computer, alignment can be conducted by block segmentation. See details.} 
}
\details{
See eRah vignette for more details. To open the vignette, execute the following code in R:
vignette("eRahManual", package="erah")

(Available also in http://metabolomicsplatform.com/applications).

For experiment of more than 1000 samples, and depending on the computer, alignment can be conducted by block segmentation. For an experiment of 1000 samples, the block.size can be set to 100, so the alignment will perform as multiple (ten) 100-samples experiments, to later align them into a single experiment.

}
\value{
The function returns an updated S4 'MetaboSet' class, where the GC-MS samples have been now aligned.
}
\references{
[1] eRah: an R package for spectral deconvolution, alignment, and metabolite identification in GC/MS-based untargeted metabolomics. Xavier Domingo-Almenara, Alexandre Perera, Maria Vinaixa, Sara Samino, Xavier Correig, Jesus Brezmes, Oscar Yanes. (2016) Article in Press.
}
\author{
Xavier Domingo-Almenara. xavier.domingo@urv.cat
}

\seealso{
\code{\link{newExp}}, \code{\link{setDecPar}}, \code{\link{deconvolveComp}}
}



