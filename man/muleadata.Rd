\name{muleaData}
\alias{muleaData}
\title{muleaData}
\usage{
muleaData()
}
\description{
Reads and downloads muleaData data.frames created from GMT files with the mulea::read_gmt() function for the mulea R package.
}
\examples{
muleaData[["Transcription_factor_TFLink_Caenorhabditis_elegans_All_EnsemblID"]]
}
\value{
Returns a data.frame that can be processed with mulea::ora, mulea::gsea, and mulea::write_gmt functions.
}
