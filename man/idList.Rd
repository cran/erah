\name{idList}
\alias{idList}
\title{Identification list}
\description{The list of identified metabolites in a given experiment}
\usage{
idList(object, id.database=mslib)
}
\arguments{
  \item{object}{A 'MetaboSet' S4 object containing the experiment data. The experiment has to be previously deconvolved, aligned and identified.}
    \item{id.database}{The mass-spectra library to be compared with the empirical spectra. By default, the MassBank - Mass Bank of North America (MoNa) database are employed (mslib object).}
}
\details{
Returns an identification table containing the names, match scores, and other variables for a given experiment.
}

\value{
	
\code{idList} returns an S3 object:
	\item{AlignID}{The unique Tag for found metabolite by eRah. Each metabolite found by eRah for a given experiment has an unique AlignID tag number.}
	\item{tmean}{The mean compound retention time.}
	\item{Name.X}{the name of the Xst/nd/rd... hit. idList return as many X (hits) as n.putative selected with \code{\link{identifyComp}}.}
	\item{FoundIn}{The number of samples in which the compound has been detected (the number of samples where the compound area is non-zero).}
	\item{MatchFactor.X}{The match factor/score of spectral similarity (spectral correlation).}
	\item{DB.Id.X}{The identification number of the library. Each metbolite in the reference library has a different DB.Id number.}
	\item{CAS.X}{the CAS number of each identified metabolite.}
}

\seealso{
\code{\link{alignList}}, \code{\link{dataList}}

}