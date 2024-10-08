% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes.R
\name{extract_model_info}
\alias{extract_model_info}
\title{Extract parameter summary from the fitted stan object}
\usage{
extract_model_info(fit, pattern = "^PS_mu|^rho|^cross|fixed")
}
\arguments{
\item{fit}{rstan fit object: fitted stan object}

\item{pattern}{string: a regex pattern corresponding to the desired parameters,
default is "^PS_mu|^rho|^cross|fixed", which are all the fixed parameters}
}
\value{
a matrix of the model output summary of the specified parameters
}
\description{
Extract parameter summary from the fitted stan object
}
