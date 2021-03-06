% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LAR.R
\name{LAR}
\alias{LAR}
\title{Estimate LAR with Single ID}
\usage{
LAR(
  data = data,
  weight_site = "no",
  weight_value = 0,
  current = as.numeric(substr(Sys.time(), 1, 4)),
  sim = 300,
  seed = 99,
  basepy = 1e+05,
  DDREF = TRUE,
  excel = FALSE,
  filename = NULL,
  ci = 0.9,
  changedata = FALSE,
  dbaseline = 0,
  dincidence = 0,
  rounddigit = 4
)
}
\arguments{
\item{data}{data frame containing 'id', 'sex', 'birth', 'dosedist', 'dose1', 'dose2', 'dose3', 'site', 'exposure_rate'. See 'Details'.}

\item{weight_site}{vector containing the name of cancer sites to give weights.}

\item{weight_value}{numeric vector containing the value between 0 and 1 which is a weight on ERR model.}

\item{current}{number of current year. default is year of the system time.}

\item{sim}{number of iteration of simulation.}

\item{seed}{number of seed.}

\item{basepy}{number of base person year}

\item{DDREF}{logical. Whether to apply the dose and dose-rate effectiveness factor.}

\item{excel}{logical. Whether to extract the result as csv file.}

\item{filename}{a string naming the file to save (.csv file).}

\item{ci}{confidence level of the confidence interval.}

\item{changedata}{logical, whether to change the data of lifetime table and incidence rate.}

\item{dbaseline}{a path or data frame of the new lifetime table.}

\item{dincidence}{a path or data frame of the new incidence rate table.}

\item{rounddigit}{the number of decimal points to print.}
}
\value{
\code{LAR}

Cancer incidence probability per 100,000 persons to radiation exposure for their lifetime after exposed  year.

\code{LBR}

Lifetime baseline risk. Cumulative baseline probability of having cancer over the maximum lifetime without radiation exposure after  exposed year.

\code{LFR}

Lifetime fractional risk. Ratio LAR/LBR.

\code{Future_LAR}

LAR after current year.

\code{BFR}

Baseline future risk. Cumulative baseline probability of having cancer over the maximum lifetime without radiation exposure after  current year.

\code{TFR}

Total future risk. Future LAR + BFR
}
\description{
\code{LAR} is used to estimate lifetime attributable radiation-related cancer risk for data with single ID.
}
\details{
The data to be put in \code{LAR} should include some prerequisite information, which includes id, sex and birth of people(or person), distribution of dose, doses of exposed radiation, sites where exposed, and exposure rate.
Also, the variable names should be written as expressed.
The maximum age in \code{LAR} is set as 100 years old. If the data contains a birth year which makes attained age(= current - birth) over 100, the result has no useful value.
For some variables, there is a fixed format. \code{sex} can have the component 'male' or 'female'.
\code{dosedist} can have the component 'fixedvalue', 'normal', 'lognormal', 'triangular', 'logtriangular', 'uniform', 'loguniform'.
\code{site} can have the component 'stomach', 'colon', 'liver', 'lung', 'breast', 'ovary', 'uterus', 'prostate', 'bladder', 'brain/cns',
'thyroid', 'remainder', 'oral', 'oesophagus', 'rectum', 'gallbladder', 'pancreas', 'kidney', 'leukemia'.
\code{exposure_rate} can have the component 'acute' or 'chronic'.
}
\examples{
data<-data.frame(id='a001', birth=1970, exposure=1980, dosedist='fixedvalue',
                dose1=10, dose2=0, dose3=0, sex='male', site='colon',
                exposure_rate='acute')

LAR(data)
}
\references{
Berrington de Gonzalez, A., Iulian Apostoaei, A., Veiga, L., Rajaraman, P., Thomas, B., Owen Hoffman, F., Gilbert, E. and Land, C. (2012). RadRAT: a radiation risk assessment tool for lifetime cancer risk projection. \emph{Journal of Radiological Protection}, \bold{32(3)}, pp.205-222.

National Research Council (NRC) and Committee to Assess Health Risks from Exposure to Low Levels of Ionizing Radiation (2005) \emph{Health Risks from Exposure to Low Levels of Ionizing Radiation: BEIR VII Phase 2} (Washington, DC: National Academy of Sciences)
}
