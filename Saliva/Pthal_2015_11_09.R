library(data.table)
library(StatMatch)
library(readxl)
library(Hmisc)
library(ggplot2)
library(reshape2)
library(mgcv)
library("colorspace", lib.loc="~/R/win-library/3.1")
library(gee)
library(foreign)
library(stats)
library(rmarkdown)

Pthal<-phwide2_2015.11.08

sala<-Saliva20150916
SPNov<-merge(sala, Pthal, by = "folio")
names(SPNov)

mod_BzP3T<-gam(gest_age_days_d~s(log(MBzP_3T),fx=FALSE),
               na.action=na.omit, data=SPNov)
plot(mod_BzP3T,scale=0, col=12, main ='Third Trimester MBzP', ylab="Difference in Gestational Age (days)")
summary(mod_BzP3T)

mod_BzP3Ta<-gam(gest_age_days_d~s(log(MBzP_3T),fx=TRUE, k=5),
               na.action=na.omit, data=SPNov)

plot(mod_BzP3Ta,scale=0, col=12, main ='Third Trimester MBzP', ylab="Difference in Gestational Age (days)")
summary(mod_BzP3Ta)

lmBzP3T<-lm(gest_age_days_d~log(MBzP_3T), data=SPNov)
plot(lmBzP3T)
summary(lmBzP3T)
MZ3<-plot(log(SPNov$MBzP_3T),SPNov$gest_age_days_d)
#######MEHHP_3T
mod_2Pht<-gam(gest_age_days_d~s(log(MBzP_3T))+s(log(MEHHP_3T)),
               na.action=na.omit, data=SPNov)
plot(mod_2Pht,scale=0, col=12, main ='Third Trimester 2', ylab="Difference in Gestational Age (days)")
summary(mod_2Pht)
lm2Pht<-lm(gest_age_days_d~log(MBzP_3T)+log(MEHHP_3T), data=SPNov)
plot(lm2Pht)
summary(lm2Pht)
#######MBZP 2T

mod_BzP2T<-gam(gest_age_days_d~s(log(MBzP_2T),fx=FALSE),
               na.action=na.omit, data=SPNov)
plot(mod_BzP2T,scale=0, col=12, main ='Second Trimester MBzP', ylab="Difference in Gestational Age (days)")
summary(mod_BzP2T)

lmBzP2T<-lm(gest_age_days_d~log(MBzP_2T), data=SPNov)
plot(lmBzP2T)
summary(lmBzP2T)
MZ2<-plot(log(SPNov$MBzP_2T),SPNov$gest_age_days_d)


names(SPNov)

#######  MECPP_2T

mod_MECPP_2T<-gam(gest_age_days_d~s(log(MECPP_2T),fx=FALSE),
               na.action=na.omit, data=SPNov)
plot(mod_MECPP_2T,scale=0, col=12, main ='Second Trimester MECCP', ylab="Difference in Gestational Age (days)")
summary(mod_MECPP_2T)

lm_MECPP_2T<-lm(gest_age_days_d~log(MECPP_2T), data=SPNov)
summary(lm_MECPP_2T)

#######  MECPP_3T
mod_MECPP_3T<-gam(gest_age_days_d~s(log(MECPP_3T),fx=FALSE),
                  na.action=na.omit, data=SPNov)
plot(mod_MECPP_3T,scale=0, col=12, main ='Third Trimester MECPP', ylab="Difference in Gestational Age (days)")
summary(mod_MECPP_3T)

lm_MECPP_3T<-lm(gest_age_days_d~log(MECPP_3T), data=SPNov)
summary(lm_MECPP_3T)

############MEHHP_2T # COMPLETELY NULL
mod_MEHHP_2T<-gam(gest_age_days_d~s(log(MEHHP_2T),fx=FALSE),
                  na.action=na.omit, data=SPNov)
plot(mod_MEHHP_2T,scale=0, col=12, main ='Second  Trimester MEHHP', ylab="Difference in Gestational Age (days)")
summary(mod_MEHHP_2T)

lm_MEHHP_2T<-lm(gest_age_days_d~log(MEHHP_2T), data=SPNov)
summary(lm_MEHHP_2T)

##############MEHHP_3T   P=0.04 Linear and negatively associated 
mod_MEHHP_3T<-gam(gest_age_days_d~s(log(MEHHP_3T),fx=FALSE),
                  na.action=na.omit, data=SPNov)
plot(mod_MEHHP_3T,scale=0, col=12, main ='Third Trimester MEHHP', ylab="Difference in Gestational Age (days)")
summary(mod_MEHHP_3T)

lm_MEHHP_3T<-lm(gest_age_days_d~log(MEHHP_3T), data=SPNov)
summary(lm_MEHHP_3T)

##############MBP_2T   COMPLETELY NULL

mod_MBP_2T<-gam(gest_age_days_d~s(log(MBP_2T),fx=FALSE),
                  na.action=na.omit, data=SPNov)
plot(mod_MBP_2T,scale=0, col=12, main ='Second Trimester MBP', ylab="Difference in Gestational Age (days)")
summary(mod_MBP_2T)

lm_MBP_2T<-lm(gest_age_days_d~log(MBP_2T), data=SPNov)
summary(lm_MBP_2T)

###############MBP_3T   N=481, P=0.2 Negative but not significant, linear

mod_MBP_3T<-gam(gest_age_days_d~s(log(MBP_3T),fx=FALSE),
                na.action=na.omit, data=SPNov)
plot(mod_MBP_3T,scale=0, col=12, main ='Third Trimester MBP', ylab="Difference in Gestational Age (days)")
summary(mod_MBP_3T)

lm_MBP_3T_2<-lm(gest_age_days_d~log(MBP_3T)+log(MBzP_3T), data=SPNov)
summary(lm_MBP_3T_2)
lm_MBP_3T<-lm(gest_age_days_d~log(MBP_3T), data=SPNov)
summary(lm_MBP_3T)

##############MEOHP_2T Null, not negative
mod_MEOHP_2T<-gam(gest_age_days_d~s(log(MEOHP_2T),fx=FALSE),
                na.action=na.omit, data=SPNov)
plot(mod_MEOHP_2T,scale=0, col=12, main ='Second Trimester MEOHP', ylab="Difference in Gestational Age (days)")
summary(mod_MEOHP_2T)

lm_MEOHP_2T<-lm(gest_age_days_d~log(MEOHP_2T), data=SPNov)
summary(lm_MEOHP_2T)


###########MEOHP_3T n=481, p=0.03, linear
mod_MEOHP_3T<-gam(gest_age_days_d~s(log(MEOHP_3T),fx=FALSE),
                  na.action=na.omit, data=SPNov)
plot(mod_MEOHP_3T,scale=0, col=12, main ='Third Trimester MEOHP', ylab="Difference in Gestational Age (days)")
summary(mod_MEOHP_3T)

lm_MEOHP_3T<-lm(gest_age_days_d~log(MEOHP_3T), data=SPNov)
summary(lm_MEOHP_3T)


############MiBP_2T NULL

mod_MiBP_2T<-gam(gest_age_days_d~s(log(MiBP_2T),fx=FALSE),
                 na.action=na.omit, data=SPNov)
plot(mod_MiBP_2T,scale=0, col=12, main ='Second Trimester MiBP', ylab="Difference in Gestational Age (days)")
summary(mod_MiBP_2T)

lm_MiBP_2T<-lm(gest_age_days_d~log(MiBP_2T), data=SPNov)
summary(lm_MiBP_2T)
############MiBP_3T COMPLETELY NULL

mod_MiBP_3T<-gam(gest_age_days_d~s(log(MiBP_3T),fx=FALSE),
                  na.action=na.omit, data=SPNov)
plot(mod_MiBP_3T,scale=0, col=12, main ='Third Trimester MiBP', ylab="Difference in Gestational Age (days)")
summary(mod_MiBP_3T)

lm_MiBP_3T<-lm(gest_age_days_d~log(MiBP_3T), data=SPNov)
summary(lm_MiBP_3T)

#########MEP_2T NEGATIVE, P=0.2, linear, n=312

mod_MEP_2T<-gam(gest_age_days_d~s(log(MEP_2T),fx=FALSE),
                 na.action=na.omit, data=SPNov)
plot(mod_MEP_2T,scale=0, col=12, main ='Second Trimester MEP', ylab="Difference in Gestational Age (days)")
summary(mod_MEP_2T)

lm_MEP_2T<-lm(gest_age_days_d~log(MEP_2T), data=SPNov)
summary(lm_MEP_2T)


###########MEP_3T Largely negative, but null 0.195


mod_MEP_3T<-gam(gest_age_days_d~s(log(MEP_3T),fx=FALSE),
                na.action=na.omit, data=SPNov)
plot(mod_MEP_3T,scale=0, col=12, main ='Third Trimester MEP', ylab="Difference in Gestational Age (days)")
summary(mod_MEP_3T)

lm_MEP_3T<-lm(gest_age_days_d~log(MEP_3T), data=SPNov)
summary(lm_MEP_3T)

############MCPP_2T 3 0's so log of MCPP is infinity, negative but null, n=309

mod_MCPP_2T<-gam(gest_age_days_d~s(log(MCPP_2T),fx=FALSE),
                na.action=na.omit, data=SPNov[SPNov$MCPP_2T!=0,])
plot(mod_MCPP_2T,scale=0, col=12, main ='Second Trimester MCPP', ylab="Difference in Gestational Age (days)")
summary(mod_MCPP_2T)

lm_MCPP_2T<-lm(gest_age_days_d~log(MCPP_2T), data=SPNov[SPNov$MCPP_2T!=0,])
summary(lm_MCPP_2T)


############MCPP_3T NULL, P=0.56 on linear model, n=475

mod_MCPP_3T<-gam(gest_age_days_d~s(log(MCPP_3T),fx=FALSE),
                 na.action=na.omit, data=SPNov[SPNov$MCPP_3T!=0,])
plot(mod_MCPP_3T,scale=0, col=12, main ='Third Trimester MCPP', ylab="Difference in Gestational Age (days)")
summary(mod_MCPP_3T)

lm_MCPP_3T<-lm(gest_age_days_d~log(MCPP_3T), data=SPNov[SPNov$MCPP_3T!=0,])
summary(lm_MCPP_3T)




###############  T.Test with preterm


###MBzP_2T P=-.08, higher in preterm
t.test(log(MBzP_2T)~preterm, data=SPNov)

###MBzP_3T P=0.03, higher in preterm 434 term, 47 preterm
t.test(log(MBzP_3T)~preterm, data=SPNov)

###MBP_2T P=0.5, hihger in preterm
t.test(log(MBP_2T)~preterm, data=SPNov)
table(SPNov$preterm, SPNov$shipmon_3T)
names(SPNov)
###MBP_3T P=0.3, higher in preterm
t.test(log(MBP_3T)~preterm, data=SPNov)



###MEOHP_2T, P=0.26, lower in preterm
t.test(log(MEOHP_2T)~preterm, data=SPNov)

###MEOHP_3T P=0.12, higher in preterm
t.test(log(MEOHP_3T)~preterm, data=SPNov)

##MiBP_2T P=0.5, lower in preterm
t.test(log(MiBP_2T)~preterm, data=SPNov)

###MiBP_3T P=0.8, higher in preterm, but really null
t.test(log(MiBP_3T)~preterm, data=SPNov)

##MEP_2T P=0.27, higher in preterm
t.test(log(MEP_2T)~preterm, data=SPNov)

###MEP_3T 
t.test(log(MEP_3T)~preterm, data=SPNov)




glm_MBzP_3T<-glm(preterm~MBzP_3T_LIQR, family="binomial",data=SPNov)
summary(glm_MBzP_3T)
exp(0.4274)
summary(log(SPNov$MBzP_3T))
table(SPNov$preterm, useNA="ifany")
SPNov$MBzP_3T_LIQR<-log(SPNov$MBzP_3T)/(IQR(log(SPNov$MBzP_3T[!is.na(SPNov$preterm+SPNov$MBzP_3T)])))
summary(SPNov$MBzP_3T_LIQR)

SPNov$meanlogdhea
plot(SPNov$meanlogdhea, log(SPNov$MBzP_3T))
plot(SPNov$meanlogdhea, log(SPNov$MECPP_2T))

with(SPNov, table(!is.na(MBzP_2T), !is.na(MBzP_3T), !is.na(preterm)))
exp(mean(log(SPNov$MBzP_3T), na.rm=TRUE))
exp(mean(log(SPNov$MBzP_2T), na.rm=TRUE))

