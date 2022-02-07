#'Average Estimated Lifetime Attributable Risk by Group
#'
#'\code{LAR_group} is used to estimate lifetime attributable radiation-related cancer risk by group.
#'
#'@param group a vector or list of vectors which distinguish each group.
#'@inheritParams LAR_batch
#'
#'@return \code{LAR_group} returns an object of multiple classes
#'  "\code{LAR_group}", "\code{LAR}". An object of class \code{LAR_group} is a
#'  list of \code{LAR} class objects which names of elements are \code{group} of each
#'  groups.
#'
#'@examples
#'## example with lifetime and incidence rate table in 2010 Korea.
#'lar1 <- LAR_group(nuclear, pid=nuclear$ID, group=nuclear$distance,
#'                  basedata = list(life2010, incid2010))
#'summary(lar1)
#'
#'lar2 <- LAR_group(nuclear, pid=nuclear$ID, group=list(nuclear$sex, nuclear$distance),
#'                  basedata = list(life2010, incid2010))
#'summary(lar2)
#'
#'
#'@references Berrington de Gonzalez, A., Iulian Apostoaei, A., Veiga, L., Rajaraman, P., Thomas, B., Owen Hoffman, F., Gilbert, E. and Land, C. (2012). RadRAT: a radiation risk assessment tool for lifetime cancer risk projection. \emph{Journal of Radiological Protection}, \bold{32(3)}, pp.205-222.
#'@references National Research Council (NRC) and Committee to Assess Health Risks from Exposure to Low Levels of Ionizing Radiation (2005) \emph{Health Risks from Exposure to Low Levels of Ionizing Radiation: BEIR VII Phase 2} (Washington, DC: National Academy of Sciences)
#'
#'@export
#'@importFrom Rcpp evalCpp
#'@importFrom stats quantile
#'@importFrom stats aggregate
#'@import dplyr
#'@useDynLib LARisk, .registration = TRUE


LAR_group <- function(data, pid, group, basedata, sim=300, seed=99, current=as.numeric(substr(Sys.Date(),1,4)), ci=0.9, weight=NULL, DDREF=TRUE, basepy=1e+05){
  set.seed(seed=seed)

  ##---Data check----------------------------####
  data <- within(data, {
    sex <- tolower(sex)
    site <- tolower(site)
    exposure_rate <- tolower(exposure_rate)
    dosedist <- tolower(dosedist)
    if(!("dose2" %in% names(data))) dose2 <- NA
    if(!("dose3" %in% names(data))) dose3 <- NA
  })
  data$pid <- pid

  ## split data
  if(is.list(group)){
    data_g <- split(data, c(list(pid), group), drop=TRUE)
  }else{
    data_g <- split(data, list(pid, group), drop=TRUE)
  }

  lapply(data_g, function(dd) check_data(dd, current=current))

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

  ## calculate grouped result
  data_g <- split(data, group, drop=TRUE)

  result_list <- list()
  for(index in 1:length(data_g)){

    ##---LAR simulation using C-------------------####
    ### input data for c
    data_c <- within(data_g[[index]], {
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

    SUMMAT <- F_SUMMAT <- matrix(0, nrow=nrow(data_g[[index]]), ncol=sim)
    if(any(data_c$site==19)){
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
      if((data_c[n,]$dosedist==1) & (data_c[n,]$dose1==0)){
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
    data_b <- lapply(split(data_c, list(data_g[[index]]$site, data_c$pid), drop=TRUE), function(dat_s) dat_s[which.min(dat_s$exposure),] )

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
    site_sum <- SUMMAT %>% as.data.frame %>% split(data_g[[index]]$site) %>% lapply(function(s_mat){apply(s_mat, MARGIN=2, sum) / length(unique(data_c$pid)) } ) %>% as.vector
    F_site_sum <- F_SUMMAT %>% as.data.frame %>% split(data_g[[index]]$site) %>% lapply(function(s_mat){ apply(s_mat, MARGIN=2, sum) / length(unique(data_c$pid)) }) %>% as.vector

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
      lar$leukemia <- colSums(leukemia_MAT) / length(unique(data_g[[index]]$pid))
      F_lar$leukemia <- colSums(F_leukemia_MAT) / length(unique(data_g[[index]]$pid))
    }else{
      lar$leukemia <- rep(0,3)
      F_lar$leukemia <- rep(0,3)
    }
    names(lar$leukemia) <- names(F_lar$leukemia) <- c("Lower", "Mean", "Upper")

    lar <- lar[site_list[site_list %in% names(lar)]]
    F_lar <- F_lar[site_list[site_list %in% names(F_lar)]]

    ### solid & total
    solid_sum <- colSums(subset(SUMMAT, data_c$site!=19), na.rm=TRUE) / length(unique(data_c$pid))
    total_sum <- colSums(SUMMAT, na.rm=TRUE) / length(unique(data_c$pid))
    F_solid_sum <- colSums(subset(F_SUMMAT, data_c$site!=19), na.rm=TRUE) / length(unique(data_c$pid))
    F_total_sum <- colSums(F_SUMMAT, na.rm=TRUE) / length(unique(data_c$pid))

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
    lbr <- sapply(split(lbr, sapply(data_b, function(tmp_dat) tmp_dat$site)), sum) / length(unique(data_c$pid))
    bfr <- sapply(split(bfr, sapply(data_b, function(tmp_dat) tmp_dat$site)), sum) / length(unique(data_c$pid))

    names(lbr) <- names(bfr) <- site_list[as.numeric(names(lbr))]

    if(!any(data_c$site==19)) lbr["leukemia"] <- bfr["leukemia"] <- 0

    lbr <- lbr[site_list[site_list %in% names(lbr)]]
    bfr <- bfr[site_list[site_list %in% names(bfr)]]

    lbr[c("solid", "total")] <- c(sum(lbr[names(lbr)!='leukemia']), sum(lbr))
    bfr[c("solid", "total")] <- c(sum(bfr[names(bfr)!='leukemia']), sum(bfr))

    ### Lifetime Fractional Risk
    lfr <- sapply(lar, function(s) s["Mean"]) / lbr
    names(lfr) <- names(lar)

    ### Total Future Risk
    tfr <- sapply(F_lar, function(s) s["Mean"]) + bfr
    names(tfr) <- names(F_lar)

    pinfo <- aggregate(data_g[[index]]$birth, by=list(data_g[[index]]$sex), function(dd) c( count=length(dd), birth=mean(dd)))
    pinfo <- data.frame(sex=pinfo$Group.1, pinfo$x)

    result_list[[index]] <- structure(list(LAR=lar, F_LAR=F_lar, LBR=lbr, BFR=bfr, LFR=lfr, TFR=tfr, current=current, ci=ci, pinfo=pinfo), class="LAR")

  } ## end grouped result

  names(result_list) <- names(data_g)

  class(result_list) <- c("LAR_group", "LAR")

  return(result_list)
}

