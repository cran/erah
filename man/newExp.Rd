\name{newExp}
\alias{newExp}
\title{New Experiment}
\description{Sets a new experiment for eRah}
\usage{
newExp(instrumental, phenotype=NULL, info=character())
}
\arguments{
  \item{instrumental}{The path where the instrumental .csv file is located.}
  \item{phenotype}{(optional) The path where the phenotypic .csv file is located.}
  \item{info}{Experiment description}
}
\details{
See eRah vignette for more details. To open the vignette, execute the following code in R:
vignette("eRahManual", package="erah")

(Available also in http://metabolomicsplatform.com/applications).

}
\value{
\code{newExp} returns an S4 object of the class 'MetaboSet'.
}
\references{
[1] eRah: an R package for spectral deconvolution, alignment, and metabolite identification in GC/MS-based untargeted metabolomics. Xavier Domingo-Almenara, Alexandre Perera, Maria Vinaixa, Sara Samino, Xavier Correig, Jesus Brezmes, Oscar Yanes. (2016) Article in Press.
}
\author{
Xavier Domingo-Almenara. xavier.domingo@urv.cat
}

\examples{

# Store all the raw data files in one different folder per class,
# and all the class-folders in one folder, which is the experiment
# folder. Then execute

# createdt(path)

# where path is the experiment folder path.
# The experiment can be now started by:

# ex <- newExp(instrumental="path/DEMO_inst.csv", 
# phenotype="path/DEMO_pheno.csv", info="DEMO Experiment")


}

\seealso{
\code{\link{createdt}}, \code{\link{setDecPar}}, \code{\link{setAlPar}} 
}

