% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smesir.R
\name{smesir}
\alias{smesir}
\title{smesir: Fitting Semiparametric Mixed-Effects SIR Models}
\usage{
smesir(
  formula,
  data,
  epi_params,
  region_names = NULL,
  prior = NULL,
  chains = 4,
  iter = 50000,
  warmup = 0,
  thin = max(floor((iter - warmup)/1000), 1),
  min_adaptation_cycles = 5,
  min_samps_per_cycle = NULL,
  tempering_ratio = 0.2,
  quiet = TRUE,
  sr_style = NULL,
  seed = NULL
)
}
\arguments{
\item{formula}{Object of class "formula"}

\item{data}{A named list containing the data with which the model will be fit. 
The list should include an entry for each term in the accompanying \code{formula} argument.}

\item{epi_params}{Epidemiologic parameters which must be specified by the user (typically obtained from side-information):
\code{region_populations} - Vector of populations of the regions studied, listed in the same order in which they are indexed in the data
\code{outbreak_times} - Vector of indices of the time interval at which the first cases are reported in each region
\code{mean_removal_time} - Average amount of time (in number or fractions of time intervals) that an individual remains infectious
\code{psi} - Vector of probabilities whose element \code{i} is the probability that the response event (case detection, death) occurs an infected individual \code{i - 1} time intervals after their infection}

\item{region_names}{Vector of names of the regions studied, listed in the same order in which they are indexed in the data}

\item{prior}{(Optional - reasonable default values are specified internally) A length 4 named list containing:
\code{ell} - lengthscale of the squared exponential kernal for the temporal random effect
\code{V0} - 3 variance hyperparameters for gaussian priors on the intercepts, coefficients, and temporal random effects;
\code{IGSR} - 3 pairs of shape and rate hyperparameters for inverse-gamma priors;
\code{expected_initial_infected} - the expected size of the infected population that appears at the beginning of the outbreak, used in an exponential prior;}
}
\description{
This function fits an SIR model with local, time-varying transmission 
rates to epidemic incidence data (deaths) observed in one or 
more regions. The transmission rate in each region is modeled
as a linear function of covariates plus a smooth temporal random 
effect drawn from a Gaussian Process distribution. The intercepts,
coefficients, and even the temporal random effects from the various
regions are assumed to be independent draws from global distributions
centered at "global" intercept, coefficient, and temporal random effect
values.
}
