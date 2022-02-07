## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(LARisk)

## ---- eval=FALSE--------------------------------------------------------------
#  LAR(data, basedata, sim=300, seed=99, current=as.numeric(substr(Sys.Date(),1,4)),
#      ci=0.9, weight=NULL, DDREF=TRUE, basepy=1e+05)

## ---- error=TRUE--------------------------------------------------------------
ex_data <- data.frame(sex = 'male', birth = 1900, exposure = 1980,
                  site = 'stomach', exposure_rate = "chronic",
                  dosedist = 'fixedvalue', dose1 = 10, dose2=NA, dose3=NA)

LAR(ex_data, basedata=list(life2010, incid2010)) ## error

## ---- eval=FALSE--------------------------------------------------------------
#  LAR(data,
#      basedata = list("the first is lifetime table", "the second is cancer incidence rate table"))

## -----------------------------------------------------------------------------
head(life2010)      ## lifetime table of the Korean in 2010.

## -----------------------------------------------------------------------------
head(incid2010)     ## cancer incidence rate table of the Korean in 2010.

## ---- eval=FALSE--------------------------------------------------------------
#  LAR(data, basedata, weight=list(stomach = 0.5))

## -----------------------------------------------------------------------------
ex_data <- data.frame(sex = 'male', birth = 1990, exposure = 2015,
                  site = 'leukemia', exposure_rate = "chronic",
                  dosedist = 'fixedvalue', dose1 = 10, dose2=NA, dose3=NA)

LAR(ex_data, basedata=list(life2010, incid2010), DDREF=TRUE)
LAR(ex_data, basedata=list(life2010, incid2010), DDREF=FALSE) ## the result are same

## ---- eval=FALSE--------------------------------------------------------------
#  LAR(data, basedata, seed=1111)    ## changing seed number, the result is also changed
#  LAR(data, basedata, sim=1000)     ## the large 'sim' offers a stable simulation result
#  LAR(data, basedata, basepy=1e+03) ## setting the baseline person-year is 1000

## ---- eval=FALSE--------------------------------------------------------------
#  LAR(data, basedata, current=2019) ## setting the current year is 2019

## ---- eval=FALSE--------------------------------------------------------------
#  LAR(data, basedata, ci=0.8) ## setting the confidence level is 0.8

## -----------------------------------------------------------------------------
nuclear1 <- nuclear[nuclear$ID=="ID01",]

print(nuclear1)

LAR(nuclear1, basedata = list(life2010, incid2010))


## -----------------------------------------------------------------------------
summary(LAR(nuclear1, basedata = list(life2010, incid2010)))

## -----------------------------------------------------------------------------
ex_batch <- LAR_batch(nuclear, pid=nuclear$ID, basedata = list(life2010, incid2010))

class(ex_batch)

class(ex_batch[[1]])

## -----------------------------------------------------------------------------
print(ex_batch, max.id=3)

## -----------------------------------------------------------------------------
summary(ex_batch, max.id=3)

## -----------------------------------------------------------------------------
ex_group1 <- LAR_group(nuclear, pid = nuclear$ID, group = nuclear$distance,
                       basedata = list(life2010, incid2010))
summary(ex_group1)

## ---- eval=FALSE--------------------------------------------------------------
#  write_LAR(x, filename)

## -----------------------------------------------------------------------------
head(organ)

## -----------------------------------------------------------------------------
organ1 <- organ[organ$ID=='ID01',]
ex_organ1 <- LAR(organ1, baseda=list(life2018, incid2018), current=2021)

ex_organ1

## -----------------------------------------------------------------------------
summary(ex_organ1)

## -----------------------------------------------------------------------------
ex_organ2 <- LAR_group(organ, pid=organ$ID, group=organ$sex,
                       basedata=list(life2018, incid2018), current=2021)

summary(ex_organ2)

## -----------------------------------------------------------------------------
ex_organ3 <- LAR_group(organ, pid=organ$ID, group=list(organ$sex, organ$occup),
                       basedata=list(life2018, incid2018), current=2021)

print(ex_organ3, max.id=3)

## -----------------------------------------------------------------------------
str(nuclear)

## ---- echo=FALSE--------------------------------------------------------------
hist(nuclear$dose1, main="Exposure dose", xlab="", breaks=100)

## ---- echo=FALSE--------------------------------------------------------------
ddd <- organ[!duplicated(organ$ID), c(1:3,11)]
knitr::kable(cbind(ddd[1:10,], ddd[11:20,]),
caption = "people in organ dataset", row.names = FALSE, align='c')

## -----------------------------------------------------------------------------
str(organ)

## ---- echo=FALSE--------------------------------------------------------------
hist(organ$dose1, main="Exposure dose", xlab="", breaks=60)

