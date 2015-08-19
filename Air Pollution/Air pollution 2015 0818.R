# Work done to figure out exposure for R01
library("sas7bdat", lib.loc="~/R/win-library/3.1")
library("foreign", lib.loc="C:/Program Files/R/R-3.1.3/library")
library("grid", lib.loc="C:/Program Files/R/R-3.1.3/library")
library("stats", lib.loc="C:/Program Files/R/R-3.1.3/library")
library("nlme", lib.loc="C:/Program Files/R/R-3.1.3/library")
library("ggplot2", lib.loc="~/R/win-library/3.1")
library("mgcv", lib.loc="C:/Program Files/R/R-3.1.3/library")

MAH<-read.sas7bdat("C:/R_Work_Dir/Air Pollution/august2015.sas7bdat")
gaplot<-plot(MAH$pmtri1, MAH$gest_age_days_d)
plot(MAH$pmtri1)
plot(MAH$gest_age_days_d)
plot(MAH$pmtri2)
MAH$preterm
gaplot2<-plot(MAH$pmtri2, MAH$gest_age_days_d)
MAH$gest_age_days_d
MAH$logmn2T<-log(MAH$re_MnM2T)
plot(MAH$logmn2T)
plot(MAH$logmn2T[MAH$gest_age_days_d<250], MAH$gest_age_days_d[MAH$gest_age_days_d<250])
MAHPT<-MAH[MAH$gest_age_days_d<259,]

hist(MAHPT$gest_age_days_d)
MAHPT$preterm

# Among Preterm infants
plotPT<-plot(MAHPT$pmtri1, MAHPT$gest_age_days_d)
lmPT<-lm(gest_age_days_d~pmtri1, data=MAHPT)
summary(lmPT)
abline(lmPT)
RPT1<-resid(lmPT)
plot(RPT1)

plotPT2<-plot(MAHPT$pmtri2, MAHPT$gest_age_days_d)
lmPT2<-lm(gest_age_days_d~pmtri2, data=MAHPT)
summary(lmPT2)
abline(lmPT2)


plotPT3<-plot(MAHPT$pmtri3, MAHPT$gest_age_days_d)
lmPT3<-lm(gest_age_days_d~pmtri3, data=MAHPT)
summary(lmPT3)
abline(lmPT3)

plotPT3<-plot(MAHPT$pmtri3, MAHPT$gest_age_days_d)
lmPT3a<-lm(gest_age_days_d~pmtri3+BMI_2T, data=MAHPT)
summary(lmPT3a)
abline(lmPT3a)

MAH266<-MAH[MAH$gest_age_days_d < 266,]

plotET<-plot(MAH266$pmtri1, MAH266$gest_age_days_d)
lmET<-lm(gest_age_days_d~pmtri1, data=MAH266)
summary(lmET)


##########SPLINE

modPTs<-gam(gest_age_days_d~s(pmtri1,fx=TRUE,k=5) + multiparous+as.factor(ED_level)+s(BMI_2T, fx=TRUE,k=3)+
               smoke_house_outside, na.action=na.omit, data=MAH)
summary(modPTs)

plot(modPTs,pch=18, scale=0, col=3, main ='natural spline', ylab="GestAge")

modPTs2<-gam(gest_age_days_d~s(pmtri2,fx=TRUE,k=5) + multiparous+as.factor(ED_level)+s(BMI_2T, fx=TRUE,k=3)+
              smoke_house_outside, na.action=na.omit, data=MAH)
summary(modPTs2)

plot(modPTs2,pch=18, scale=0, col=3, main ='natural spline', ylab="GestAge")

modPTs3<-gam(gest_age_days_d~s(pmtri3,fx=TRUE,k=5) + multiparous+as.factor(ED_level)+s(BMI_2T, fx=TRUE,k=3)+
               smoke_house_outside, na.action=na.omit, data=MAH)
summary(modPTs3)

plot(modPTs3,pch=18, scale=0, col=3, main ='natural spline', ylab="GestAge")

modPTs3<-gam(gest_age_days_d~s(pmtri3,fx=TRUE,k=4) + multiparous+as.factor(ED_level)+s(BMI_2T, fx=TRUE,k=3)+
               smoke_house_outside, na.action=na.omit, data=MAH)
summary(modPTs3)

plot(modPTs3,pch=18, scale=0, col=3, main ='natural spline', ylab="GestAge")

modPTs3<-gam(gest_age_days_d~s(pmtri3,fx=TRUE,k=3) + edad+multiparous+as.factor(ED_level)+BMI_2T+
               smoke_house_outside, na.action=na.omit, data=MAH)
summary(modPTs3)

plot(modPTs3,pch=18, scale=0, col=3, main ='natural spline', ylab="GestAge")

reT3<-residuals(modPTs3)

# Lets look at T1 and T2 w/ cervix 

#PTGER2
names(MAH) 
cerv1<-plot(MAH$pmtri1, MAH$Cervix_PTGER2_Mean)
modcerv1<-gam(Cervix_PTGER2_Mean~s(pmtri1,fx=TRUE,k=3) +PAP+ edad+multiparous+as.factor(ED_level)+BMI_2T+
               smoke_house_outside, na.action=na.omit, data=MAH)
summary(modcerv1)

plot(modcerv1,pch=18, scale=0, col=2, main ='natural spline', ylab="PTGER2")


modcerv2<-gam(Cervix_PTGER2_Mean~s(pmtri2,fx=TRUE,k=3) + edad+multiparous+as.factor(ED_level)+BMI_2T+
                smoke_house_outside, na.action=na.omit, data=MAH)
summary(modcerv2)

plot(modcerv2,pch=18, scale=0, col=2, main ='natural spline', ylab="PTGER2")

#LINE1

cervLINE<-plot(MAH$pmtri1, MAH$LINEmethAve)
modcervLINE<-gam(LINEmethAve~s(pmtri1,fx=TRUE,k=3) +PAP+ edad+multiparous+as.factor(ED_level)+BMI_2T+
                smoke_house_outside, na.action=na.omit, data=MAH)
summary(modcervLINE)

plot(modcerv1,pch=18, scale=0, col=2, main ='natural spline', ylab="LINE1")


modcervLINE2<-gam(LINEmethAve~s(pmtri2,fx=TRUE,k=3) + edad+multiparous+as.factor(ED_level)+BMI_2T+
                smoke_house_outside, na.action=na.omit, data=MAH)
summary(modcervLINE2)

plot(modcervLINE2,pch=18, scale=0, col=2, main ='natural spline', ylab="LINEmethAve")

#MiRNAs of interest
?merge
mir<-mir_cervix_sorted

MAH$folio
mir$folio<-mir$Folio
MAHmir<-merge(MAH, mir, by="folio")
MAHmir$folio

names(MAHmir)
MAHmir$mir21_5p<-MAHmir$hsa_miR_21_5p
MAHmir$mir142_3p<-MAHmir$hsa_miR_142_3p
MAHmir$mir29b<-MAHmir$hsa_miR_29b_3p
MAHmir$mir30e_5p<-MAHmir$hsa_miR_30e_5p
MAHmir$mir148b<-MAHmir$hsa_miR_148b_3p
MAHmir$mir223<-MAHmir$hsa_miR_223_3p
names(MAHmir)
#mir 21
plot(MAHmir$mir21_5p, MAHmir$pmtri1)
lmMir21<-lm(mir21_5p~pmtri1, data=MAHmir)
summary(lmMir21)
plot(lmMir21)
abline(lmMir21)

plot(MAHmir$pmtri2,MAHmir$mir21_5p)
lmMir21T2<-lm(mir21_5p~pmtri2, data=MAHmir)
summary(lmMir21T2)
plot(lmMir21T2)
#mir 142

plot(MAHmir$pmtri1,MAHmir$mir142_3p)
lmMir142<-lm(mir142_3p~pmtri1, data=MAHmir)
summary(lmMir142)
plot(lmMir142)

R142<-plot(MAHmir$pmtri2,MAHmir$mir142_3p)
lmMir142T2<-lm(mir142_3p~pmtri2, data=MAHmir)
summary(lmMir142T2)
abline(lmMir142T2)
plot(lmMir142T2)

#mir223
plot(MAHmir$pmtri1,MAHmir$mir223)
lmMir223<-lm(mir223~pmtri1, data=MAHmir)
summary(lmMir223)
abline(lmMir223)

plot(MAHmir$pmtri2,MAHmir$mir223)
lmMir223T2<-lm(mir223~pmtri2, data=MAHmir)
summary(lmMir223T2)
abline(lmMir223T2)

#mir30e_5p

plot(MAHmir$pmtri1,MAHmir$mir30e_5p)
lmMir30e<-lm(mir30e_5p~pmtri1, data=MAHmir)
summary(lmMir30e)
abline(lmMir30e)

plot(MAHmir$pmtri2,MAHmir$mir30e_5p)
lmMir30eT2<-lm(mir30e_5p~pmtri2, data=MAHmir)
summary(lmMir30eT2)
abline(lmMir30eT2)



#mir148b

plot(MAHmir$pmtri1,MAHmir$mir148b)
lmMir148<-lm(mir148b~pmtri1, data=MAHmir)
summary(lmMir148)
abline(lmMir148)

plot(MAHmir$pmtri2,MAHmir$mir148b)
lmMir148T2<-lm(mir148b~pmtri2, data=MAHmir)
summary(lmMir148T2)
abline(lmMir148T2)

#mir29b


plot(MAHmir$pmtri1,MAHmir$mir29b)
lmMir29b<-lm(mir29b~pmtri1, data=MAHmir)
summary(lmMir29b)
abline(lmMir29b)

plot(MAHmir$pmtri2,MAHmir$mir29b)
lmMir29bT2<-lm(mir29b~pmtri2, data=MAHmir)
summary(lmMir29bT2)
abline(lmMir29bT2)