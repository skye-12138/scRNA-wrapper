% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_own_reference.R
\name{ref.make}
\alias{ref.make}
\title{Transform matrix to a standared reference format for next analysis}
\usage{
ref.make(mat, ref.count = TRUE)
}
\arguments{
\item{mat}{a matrix of gene expression value , colnames are classified variables , rownames are genes}

\item{ref.count}{logic value, if the expression values in mat are count data, FALSE means the mat has been lognormalized.}
}
\value{
a standard reference data for next analysis
}
\description{
Transform matrix to a standared reference format for next analysis
}
\examples{

}
