% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RI.R
\name{computeRIerror}
\alias{computeRIerror}
\title{computeRIerror}
\usage{
computeRIerror(
  Experiment,
  id.database = mslib,
  reference.list,
  ri.error.type = c("relative", "absolute"),
  plot.results = TRUE
)
}
\arguments{
\item{Experiment}{S4 object with experiment Data, Metadata and Results. Results of experiment are used to extract RT and Compound DB Id.}

\item{id.database}{Name of the preloaded database, in this case the regular db used by erah mslib}

\item{reference.list}{List with the compounds and their attributes (AlignId...)}

\item{ri.error.type}{Specify wether absolute or relative RI error is to be computed.}

\item{plot.results}{Shows the RI/RT graphic (True by default)}
}
\description{
This function uses RI of mslib database and RT of the identified compounds to discrimine proper compound identification.
}
\details{
See eRah vignette for more details. To open the vignette, execute the following code in R:
vignette("eRahManual", package="erah")
}
\examples{
\dontrun{
ex <- computeRIerror(
  ex, 
  mslib, 
  reference.list=list(AlignID = c(45,67,92,120)), 
  ri.error.type = "relative"
)
}
}
\references{
[1] Xavier Domingo-Almenara, et al., eRah: A Computational Tool Integrating Spectral Deconvolution and Alignment with Quantification and Identification of Metabolites in GC-MS-Based Metabolomics. Analytical Chemistry (2016). DOI: 10.1021/acs.analchem.6b02927
}
\seealso{
\code{\link{showRTRICurve}}
}
\author{
Xavier Domingo-Almenara. xavier.domingo@urv.cat
}
