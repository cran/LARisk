#'Check for LAR
#'
#'@rdname check.LAR
#'@export
#'@keywords internal
check_data <- function(data, current){
  ## component check
  if(!all(data$sex %in% c("male", "female")))
    stop("'sex' has an invalid component")
  if(!all(data$site %in% c('stomach', 'colon', 'liver', 'lung', 'breast', 'ovary', 'uterus', 'prostate', 'bladder', 'brain/cns', 'thyroid', 'remainder', 'oral', 'oesophagus', 'rectum', 'gallbladder', 'pancreas', 'kidney', 'leukemia')))
    stop("'site' has an invalid component")
  if(!all(data$exposure_rate %in% c("chronic", "acute")))
    stop("'exposure_rate' has an invalid component")
  if(!all(data$dosedist %in% c('fixedvalue', 'lognormal', 'normal', 'triangular', 'logtriangular', 'uniform', 'loguniform')))
    stop("'dosedist' has an invalid component")

  if(!is.numeric(data$birth))
    stop("'birth' must be a numeric")
  if(!is.numeric(data$exposure))
    stop("'exposure' must be a numeric")
  if(!is.numeric(c(data$dose1, data$dose2, data$dose3)))
    stop("'dose' must be a numeric")

  ## checking one person
  if(length(unique(data$sex))!=1)
    stop("'sex' of a single person must be the same value")
  if(length(unique(data$birth))!=1)
    stop("'birth' of a single person must be the same value")

  ## logical check
  if(any((data$sex=="male")&(data$site %in% c('uterus', 'breast', 'ovary'))))
    stop("'uterus', 'breast' and 'ovary' cannot be chosen for male")
  if(any((data$sex=="female")&(data$site=='prostate')))
    stop("'prostate' cannot be chosen for 'male'")

  if(any((current - data$birth)>100))
    stop("Age is not allowed to be greater than 100 years.")
  if(any(data$birth>data$exposure))
    stop("'birth' cannot exceed 'exposure'")

  if(any((data$dosedist == "fixedvalue") & is.na(data$dose1) ) )
    stop("'fixedvalue' require 'dose1'")
  if(any((data$dosedist == "lognormal") & is.na(data$dose1) & is.na(data$dose2) ) )
    stop("'lognormal' require 'dose1' and 'dose2'")
  if(any((data$dosedist == "normal") & is.na(data$dose1) & is.na(data$dose2) ) )
    stop("'normal' require 'dose1' and 'dose2'")
  if(any((data$dosedist == "triangular") & is.na(data$dose1) & is.na(data$dose2) & is.na(data$dose3) ) )
    stop("'triangular' require 'dose1', 'dose2' and 'dose3'")
  if(any((data$dosedist == "logtriangular") & is.na(data$dose1) & is.na(data$dose2) & is.na(data$dose3) ) )
    stop("'logtriangular' require 'dose1', 'dose2' and 'dose3'")
  if(any((data$dosedist == "uniform") & is.na(data$dose1) & is.na(data$dose2) ) )
    stop("'uniform' require 'dose1' and 'dose2'")
  if(any((data$dosedist == "loguniform") & is.na(data$dose1) & is.na(data$dose2) ) )
    stop("'loguniform' require 'dose1' and 'dose2'")
}

#'@rdname check.LAR
#'@import dplyr
#'@export
#'@keywords internal
check_basedata <- function(basedata){
  ## check data
  if(!all(c("Age", "Prob_d_m", "Prob_d_f") %in% names(basedata[[1]])))
    stop("Lifetime table must include 'Age', 'Prob_d_m' and 'Prob_d_f'")
  if(!all(c("Site", "Age", "Rate_m", "Rate_f") %in% names(basedata[[2]])))
    stop("Cancer incidence rate table must include 'Site', 'Age', 'Rate_m' and 'Rate_f'")

  ## sort
  lifeTable <- basedata[[1]] %>% select("Age", "Prob_d_m", "Prob_d_f") %>% arrange(basedata[[1]]$Age)
  incidTable <- basedata[[2]] %>%
    mutate(Site=recode(basedata[[2]]$Site,
                       'stomach'    =1,  'colon'    =2,  'liver'   =3,  'lung'      =4,  'breast'   =5,
                       'ovary'      =6,  'uterus'   =7,  'prostate'=8,  'bladder'   =9,  'brain/cns'=10,
                       'thyroid'    =11, 'remainder'=12, 'oral'    =13, 'oesophagus'=14, 'rectum'   =15,
                       'gallbladder'=16, 'pancreas' =17, 'kidney'  =18, 'leukemia'  =19)) %>%
    select("Site", "Age", "Rate_m", "Rate_f")
  incidTable <- incidTable %>% arrange(incidTable$Site, incidTable$Age)

  if(any(lifeTable$Age != 0:100))
    stop("The range of 'Age' are from 0 to 100 by one by one")

  if(any((incidTable$Site != rep(1:19, each=101)), any(incidTable$Age != rep(0:100, 19))))
    stop("The range of 'Age' are from 0 to 100 for each 'Site' by one by one")

  return(list(lifeTable=lifeTable, incidTable=incidTable))
}
