#'Estimate Lifetime Attributable Risk for several people
#'
#'\code{LAR_batch} is used to estimate lifetime attributable radiation-related cancer risk for data with several people.
#'
#'@param pid a vector which distinguish each person.
#'@inheritParams LAR
#'
#'@return \code{LAR_batch} returns an object of multiple classes
#'  "\code{LAR_batch}", "\code{LAR}". An object of class \code{LAR_batch} is a
#'  list of \code{LAR} class objects which names of elements are \code{ID} of each
#'  person.
#'
#'@examples
#'## example with lifetime and incidence rate table in 2010 Korea.
#'lar1 <- LAR_batch(nuclear, pid=nuclear$ID, basedata = list(life2010, incid2010))
#'summary(lar1)
#'
#'@seealso \code{\link{LAR}}, \code{\link{LAR_group}}
#'
#'@references Berrington de Gonzalez, A., Iulian Apostoaei, A., Veiga, L., Rajaraman, P., Thomas, B., Owen Hoffman, F., Gilbert, E. and Land, C. (2012). RadRAT: a radiation risk assessment tool for lifetime cancer risk projection. \emph{Journal of Radiological Protection}, \bold{32(3)}, pp.205-222.
#'@references National Research Council (NRC) and Committee to Assess Health Risks from Exposure to Low Levels of Ionizing Radiation (2005) \emph{Health Risks from Exposure to Low Levels of Ionizing Radiation: BEIR VII Phase 2} (Washington, DC: National Academy of Sciences)
#'
#'@export
#'@importFrom Rcpp evalCpp
#'@importFrom stats quantile
#'@import dplyr
#'@useDynLib LARisk, .registration = TRUE

LAR_batch <- function(data, pid, basedata, sim=300, seed=99, current=as.numeric(substr(Sys.Date(),1,4)), ci=0.9, weight=NULL,DDREF=TRUE, basepy=1e+05){
  set.seed(seed = seed)

  ##---Data check----------------------------####
  data <- within(data, {
    sex <- tolower(sex)
    site <- tolower(site)
    exposure_rate <- tolower(exposure_rate)
    dosedist <- tolower(dosedist)
    if(!("dose2" %in% names(data))) dose2 <- NA
    if(!("dose3" %in% names(data))) dose3 <- NA
  })

  result <- lapply(split(data,pid), function(ss) LAR(data=ss, basedata=basedata, sim=sim, seed=seed, current=current, ci=ci, weight=weight, DDREF=DDREF, basepy=basepy))

  class(result) <- c("LAR_batch", "LAR")

  return(result)
}

