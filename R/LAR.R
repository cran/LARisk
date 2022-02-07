#'Estimate Lifetime Attributable Risk for one person
#'
#'\code{LAR} is used to estimate lifetime attributable radiation-related cancer risk for data with one person.
#'
#'@param data data frame containing demographic information and exposure information. See 'Details'.
#'@param basedata a list of the data of lifetime table and incidence rate table.
#'  The first element is lifetime table and the second is incidence rate table.
#'@param sim number of iteration of simulation.
#'@param seed a random seed number.
#'@param current a current year. default is year of the system time.
#'@param ci confidence level of the confidence interval.
#'@param weight a list containing the value between 0 and 1 which is a weight on ERR model. See 'Details'.
#'@param DDREF logical. Whether to apply the dose and dose-rate effectiveness factor.
#'@param basepy number of base person-years
#'
#'@details The maximum age in \code{LAR} is set as 100. If the data contains
#'  \code{birth} which makes attained age (=\code{current} - \code{birth})
#'  exceed 100, the result has no useful value.
#'
#'\code{data} should include information which includes gender, year of birth,
#'year of exposure, sites where exposed, exposure rate, distribution of dose and
#'dose parameters of exosed radiation. The name of each variables must be
#'\code{sex}, \code{birth}, \code{exposure}, \code{site}, \code{exposure_rate},
#'\code{dosedist}, \code{dose1}, \code{dose2}, \code{dose3}.
#'
#'For some variables, there is a fixed format. \code{sex} can have the component 'male' or 'female'.
#'\code{site} can have the component 'stomach', 'colon', 'liver', 'lung', 'breast', 'ovary', 'uterus', 'prostate', 'bladder', 'brain/cns',
#''thyroid', 'remainder', 'oral', 'oesophagus', 'rectum', 'gallbladder', 'pancreas', 'kidney', 'leukemia'.
#'\code{exposure_rate} can have the component 'acute' or 'chronic'.
#'\code{dosedist} can have the component 'fixedvalue', 'lognormal', 'normal', 'triangular', 'logtriangular', 'uniform', 'loguniform'.
#'
#'\code{dose1}, \code{dose2}, \code{dose3} are parameters of dose distribution. The parameters for each distribution are that:
#'\describe{
#'\item{fixedvalue}{dose value (dose1)}
#'\item{lognormal}{median (dose1), geometric standard deviation (dose2)}
#'\item{normal}{mean (dose1), standard deviation (dose2)}
#'\item{triangular or logtriangular}{minimum (dose1), mode (dose2), maximum (dose3)}
#'\item{uniform or loguniform}{minimum (dose1), maximum (dose2)}
#'}
#'
#'\code{weight}
#'
#'
#'@return \code{LAR} returns an object of "\code{LAR}" class.
#'
#'An object of class "\code{LAR}" is a list containing the following components:
#'\describe{
#'\item{\code{LAR}}{Lifetime attributable risk (LAR) from the time of exposure to the end of the expected lifetime.}
#'\item{\code{F_LAR}}{Future attributable risk from current to the expected lifetime.}
#'\item{\code{LBR}}{Lifetime baseline risk.}
#'\item{\code{BFR}}{Baseline future risk.}
#'\item{\code{LFR}}{Lifetime fractional risk.}
#'\item{\code{TFR}}{Total future risk.}
#'\item{\code{current}}{Current year.}
#'\item{\code{ci}}{Confidence level.}
#'\item{\code{pinfo}}{Information of the person.}
#'}
#'
#'@examples
#'## example with lifetime and incidence rate table in 2010 Korea.
#'organ2 <- split(organ, organ$ID)[[1]]   ## data of one person.
#'
#'## defualt
#'lar1 <- LAR(organ2, basedata = list(life2010, incid2010))
#'summary(lar1)
#'
#'## change the weight for ERR and EAR models
#'weight_list <- list("rectum" = 0.5)
#'lar2 <- LAR(organ2, basedata = list(life2010, incid2010), weight = weight_list)
#'summary(lar2)
#'
#'## change the DDREF option (DDREF=FALSE)
#'lar3 <- LAR(organ2, basedata = list(life2010, incid2010), DDREF = FALSE)
#'summary(lar3)
#'
#'
#'@seealso \code{\link{LAR_batch}}, \code{\link{LAR_group}}
#'
#'
#'@references Berrington de Gonzalez, A., Iulian Apostoaei, A., Veiga, L.,
#'  Rajaraman, P., Thomas, B., Owen Hoffman, F., Gilbert, E. and Land, C.
#'  (2012). RadRAT: a radiation risk assessment tool for lifetime cancer risk
#'  projection. \emph{Journal of Radiological Protection}, \bold{32(3)},
#'  pp.205-222.
#'@references National Research Council (NRC) and Committee to Assess Health
#'  Risks from Exposure to Low Levels of Ionizing Radiation (2005) \emph{Health
#'  Risks from Exposure to Low Levels of Ionizing Radiation: BEIR VII Phase 2}
#'  (Washington, DC: National Academy of Sciences)
#'
#'
#'@export
#'@importFrom Rcpp evalCpp
#'@importFrom stats quantile
#'@import dplyr
#'@useDynLib LARisk, .registration = TRUE

LAR <- function(data, basedata, sim=300, seed=99, current=as.numeric(substr(Sys.Date(),1,4)), ci=0.9, weight=NULL, DDREF=TRUE, basepy=1e+05){
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

  check_data(data, current)

  ##---Setting Global Variables----------------------####
  ### life & incidence table
  basedata <- check_basedata(basedata)
  lifeTable <- basedata$lifeTable
  incidTable <- basedata$incidTable

  ### weight for ERR models
  weight_value <- c(
    "stomach"    =0.7, "colon"    =0.7, "liver"   =0.7, "lung"      =0.3, "breast"   =0.0,
    "ovary"      =0.7, "uterus"   =0.7, "prostate"=0.7, "bladder"   =0.7, "brain/cns"=1.0,
    "thyroid"    =1.0, "remainder"=0.7, "oral"    =0.7, "oesophagus"=0.7, "rectum"   =0.7,
    "gallbladder"=1.0, "pancreas" =0.7, "kidney"  =0.7, "leukemia"  =0.7
  )
  if(!is.null(weight)) weight_value[names(weight)] <- unlist(weight)


  ##---LAR simulation using C-------------------####
  ### input data for c
  data_c <- within(data, {
    sex <- ifelse(sex=='male', 1, 2)
    birth <- birth
    exposure <- exposure
    site <- recode(site,
                   'stomach'    =1,  'colon'    =2,  'liver'   =3,  'lung'      =4,  'breast'   =5,
                   'ovary'      =6,  'uterus'   =7,  'prostate'=8,  'bladder'   =9,  'brain/cns'=10,
                   'thyroid'    =11, 'remainder'=12, 'oral'    =13, 'oesophagus'=14, 'rectum'   =15,
                   'gallbladder'=16, 'pancreas' =17, 'kidney'  =18, 'leukemia'  =19)
    exposure_rate <- ifelse(exposure_rate=='acute', 1, 2)
    dosedist <- recode(dosedist,
                       'fixedvalue'   =1, 'lognormal'=2, 'normal'    =3, 'triangular'=4,
                       'logtriangular'=5, 'uniform'  =6, 'loguniform'=7)
    age <- current - birth
    exposeage <- exposure - birth
    weight_err <- weight_value[site]
  })

  SUMMAT <- F_SUMMAT <- matrix(0, nrow=nrow(data), ncol=sim)
  if(any(data_c$site==19)){ ## exist leukemia at least one
    leukemia_MAT <- F_leukemia_MAT <- matrix(0, nrow=sum(data_c$site==19), ncol=3)
    n_leukemia <- 0
  }
  for(n in 1:nrow(data_c)){
    doseinfo <- with(data_c[n,], switch(dosedist,
                                        '1' = dose1,
                                        '2' = c(dose1, dose2),
                                        '3' = c(dose1, dose2),
                                        '4' = c(dose1, dose2, dose3),
                                        '5' = c(dose1, dose2, dose3),
                                        '6' = c(dose1, dose2),
                                        '6' = c(dose1, dose2) ))

    ## calculate results from C
    if((data_c[n,]$dosedist==1)&(data[n,]$dose1==0)){
      SUMMAT[n,] <- F_SUMMAT[n,] <-  numeric(sim)

      if(data_c$site[n]==19){
        n_leukemia <- n_leukemia + 1
        leukemia_MAT[n_leukemia,] <- F_leukemia_MAT[n_leukemia,] <- numeric(3)
      }
    }else{
      result <- .C("larft",
                   as.integer(data_c$sex[n]), as.integer(data_c$age[n]), as.integer(data_c$exposeage[n]), as.integer(data_c$site[n]),
                   as.integer(data_c$exposure_rate[n]), as.integer(data_c$dosedist[n]), as.double(doseinfo), as.double(data_c$weight_err[n]),
                   as.integer(DDREF), as.double(lifeTable$Prob_d_m), as.double(lifeTable$Prob_d_f),
                   as.double(incidTable$Rate_m), as.double(incidTable$Rate_f), as.integer(sim), as.double(ci),
                   LAR_total = as.double(vector(mode="numeric", length=sim)), F_LAR_total = as.double(vector(mode="numeric", length=sim)),
                   leukemia_result = as.double(vector(mode="numeric", length=6))
      )

      SUMMAT[n,] <- result$LAR_total*basepy/(10^5)
      F_SUMMAT[n,] <- result$F_LAR_total*basepy/(10^5)

      if(data_c$site[n]==19){
        n_leukemia <- n_leukemia + 1
        leukemia_MAT[n_leukemia,] <- result$leukemia_result[1:3]*basepy/(10^5)
        F_leukemia_MAT[n_leukemia,] <- result$leukemia_result[4:6]*basepy/(10^5)
      }
    }
  }## end LAR simulation

  ##---Baseline risk calculation using C--------####
  data_b <- lapply(split(data_c, data$site), function(dat_s) dat_s[which.min(dat_s$exposure),] )

  lbr <- unlist(lapply(data_b, function(dat){
    result <- .C('brft',
                 as.integer(dat$sex), as.integer(dat$site), as.integer(dat$exposeage),
                 as.double(lifeTable$Prob_d_m), as.double(lifeTable$Prob_d_f), as.double(incidTable$Rate_m), as.double(incidTable$Rate_f),
                 BR = as.double(vector(mode="numeric", length=1)))$BR
    return(result*basepy/(10^5))
  }))

  bfr <- unlist(lapply(data_b, function(dat){
    result <- .C('brft',
                 as.integer(dat$sex), as.integer(dat$site), as.integer(dat$age),
                 as.double(lifeTable$Prob_d_m), as.double(lifeTable$Prob_d_f), as.double(incidTable$Rate_m), as.double(incidTable$Rate_f),
                 BR = as.double(vector(mode="numeric", length=1)))$BR
    return(result*basepy/(10^5))
  }))

  ##---Cleanup Result---------------------------####
  site_list <- c('stomach', 'colon', 'liver', 'lung', 'breast', 'ovary', 'uterus', 'prostate', 'bladder', 'brain/cns', 'thyroid', 'remainder', 'oral', 'oesophagus', 'rectum', 'gallbladder', 'pancreas', 'kidney', 'leukemia')

  ### site
  site_sum <- SUMMAT %>% as.data.frame %>% split(data$site) %>% lapply(function(s_mat) apply(s_mat, MARGIN=2, sum)) %>% as.vector
  F_site_sum <- F_SUMMAT %>% as.data.frame %>% split(data$site) %>% lapply(function(s_mat) apply(s_mat, MARGIN=2, sum)) %>% as.vector

  lar <- lapply(site_sum, function(ss){
    result <- c(quantile(ss, (1-ci)/2, na.rm=TRUE), mean(ss, na.rm=TRUE), quantile(ss, (1+ci)/2, na.rm=TRUE))
    names(result) <- c("Lower", "Mean", "Upper")
    return(result)
  } )

  F_lar <- lapply(F_site_sum, function(ss){
    result <- c(quantile(ss, (1-ci)/2, na.rm=TRUE), mean(ss, na.rm=TRUE), quantile(ss, (1+ci)/2, na.rm=TRUE))
    names(result) <- c("Lower", "Mean", "Upper")
    return(result)
  } )

  if(any(data_c$site==19)){
    lar$leukemia <- colSums(leukemia_MAT)
    F_lar$leukemia <- colSums(F_leukemia_MAT)
  }else{
    lar$leukemia <- rep(0,3)
    F_lar$leukemia <- rep(0,3)
  }
  names(lar$leukemia) <- names(F_lar$leukemia) <- c("Lower", "Mean", "Upper")

  lar <- lar[site_list[site_list %in% names(lar)]]
  F_lar <- F_lar[site_list[site_list %in% names(F_lar)]]

  ### solid & total
  solid_sum <- colSums(subset(SUMMAT, data$site!='leukemia'), na.rm=TRUE)
  total_sum <- colSums(SUMMAT, na.rm=TRUE)
  F_solid_sum <- colSums(subset(F_SUMMAT, data$site!='leukemia'), na.rm=TRUE)
  F_total_sum <- colSums(F_SUMMAT, na.rm=TRUE)

  lar$solid <- c(quantile(solid_sum, (1-ci)/2, na.rm=TRUE), mean(solid_sum, na.rm=TRUE), quantile(solid_sum, (1+ci)/2, na.rm=TRUE))
  F_lar$solid <- c(quantile(F_solid_sum, (1-ci)/2, na.rm=TRUE), mean(F_solid_sum, na.rm=TRUE), quantile(F_solid_sum, (1+ci)/2, na.rm=TRUE))

  if(all(data_c$site==19)){
    lar$total <- lar$leukemia
    F_lar$total <- F_lar$leukemia
  }else{
    lar$total <- c(quantile(total_sum, (1-ci)/2, na.rm=TRUE), mean(total_sum, na.rm=TRUE), quantile(total_sum, (1+ci)/2, na.rm=TRUE))
    F_lar$total <- c(quantile(F_total_sum, (1-ci)/2, na.rm=TRUE), mean(F_total_sum, na.rm=TRUE), quantile(F_total_sum, (1+ci)/2, na.rm=TRUE))
  }

  names(lar$solid) <- names(lar$total) <- names(F_lar$solid) <- names(F_lar$total) <- c("Lower", "Mean", "Upper")

  ### Baseline Risk
  lbr <- lbr[site_list[site_list %in% names(lbr)]]
  bfr <- bfr[site_list[site_list %in% names(bfr)]]

  if(!any(data_c$site==19)) lbr["leukemia"] <- bfr["leukemia"] <- 0

  lbr[c("solid", "total")] <- c(sum(lbr[names(lbr)!='leukemia']), sum(lbr))
  bfr[c("solid", "total")] <- c(sum(bfr[names(bfr)!='leukemia']), sum(bfr))

  ### Lifetime Fractional Risk
  lfr <- sapply(lar, function(s) s["Mean"]) / lbr
  names(lfr) <- names(lar)

  ### Total Future Risk
  tfr <- sapply(F_lar, function(s) s["Mean"]) + bfr
  names(tfr) <- names(F_lar)

  pinfo <- data[1,c('sex', 'birth')]
  return(structure(list(LAR=lar, F_LAR=F_lar, LBR=lbr, BFR=bfr, LFR=lfr, TFR=tfr, current=current, ci=ci, pinfo=pinfo), class="LAR"))
}
