library("sas7bdat", lib.loc="~/R/win-library/3.1")
library("foreign", lib.loc="C:/Program Files/R/R-3.1.3/library")
library("grid", lib.loc="C:/Program Files/R/R-3.1.3/library")
library("stats", lib.loc="C:/Program Files/R/R-3.1.3/library")
library("nlme", lib.loc="C:/Program Files/R/R-3.1.3/library")
library("ggplot2", lib.loc="~/R/win-library/3.1")
library("mgcv", lib.loc="C:/Program Files/R/R-3.1.3/library")


SalivaPart<-read.sas7bdat("C:/R_Work_Dir/Saliva/salivapart.sas7bdat")
bob6<-read.sas7bdat("C:/R_Work_Dir/Saliva/bob5.sas7bdat")
names(bob6)

mergex<-read.sas7bdat("C:/R_Work_Dir/Saliva/merge.sas7bdat")
names(mergex)
#modeling w/ ahrr

SP2<-SalivaPart[SalivaPart$ahrr_mean>20,]
names(SP2)
plot(SP2$meanlogprog,SP2$ahrr_mean)


modahrr<-lm(ahrr_mean~meanlogprog, data=SP2)
plot(modahrr)
summary(modahrr)

RESAHRR<-resid(modahrr)
#AHRR and saliva prog has 194 in model

modahrr2<-lm(ahrr_mean~meanlogprog+BMI_2T+smoke_house_outside+exsmoker+
               Fenton_Z_score+edad+gest_age_days_d+multiparous+ed_gt_12+ed_LT_12, data=SP2)
summary(modahrr2)

modahrr2est<-lm(ahrr_mean~meanlogest+BMI_2T+smoke_house_outside+exsmoker+
               Fenton_Z_score+edad+gest_age_days_d+multiparous+ed_gt_12+ed_LT_12, data=SP2)
summary(modahrr2est)


#model BWT

modBWT<-lm(Fenton_Z_score~meanlogprog+edad+BMI_2T, data=SalivaPart)
summary(modBWT)

modBWTEST<-lm(Fenton_Z_score~meanlogest+edad+BMI_2T, data=SalivaPart)
summary(modBWTEST)

modBWTtEST<-lm(Fenton_Z_score~meanlogtest+edad+BMI_2T, data=SalivaPart)
summary(modBWTtEST)

modBWTmel<-lm(Fenton_Z_score~meanlogmel+edad+BMI_2T, data=SalivaPart)
summary(modBWTmel)

modBWTcort<-lm(Fenton_Z_score~meanlogcort+edad+BMI_2T, data=SalivaPart)
summary(modBWTcort)


SalivaPart$newlogcortAM<-(SalivaPart$newlogcort1+SalivaPart$newlogcort2)/2
mean(SalivaPart$newlogcortAM, na.rm=TRUE)
SalivaPart$cortslope<-(SalivaPart$newlogcortAM-SalivaPart$newlogcort3)

#No effect of cortisol on BWT
modBWTcort<-lm(Fenton_Z_score~cortslope+edad+BMI_2T, data=SalivaPart)
summary(modBWTcort)

modBWTdhea<-lm(Fenton_Z_score~meanlogdhea+edad+BMI_2T, data=SalivaPart)
summary(modBWTdhea)


# GEST AGE

modGAprog<-lm(gest_age_days_d~meanlogprog, data=SalivaPart) #n=219
summary(modGAprog)
RRGA<-resid(modGAprog)

modPreterm<-lm(meanlogprog~preterm,data=SalivaPart)
summary(modPreterm)

?stats
library(help="stats")
modprog<-lm(meanlogprog~male+BMI_2T+ed_LT_12+ed_gt_12+multiparous+
              smoke_house_outside, data=SalivaPart)
summary(modprog)

modest<-lm(meanlogest~male+BMI_2T+ed_LT_12+ed_gt_12+multiparous+
              smoke_house_outside, data=SalivaPart)
summary(modest)

modtest<-lm(meanlogtest~male+BMI_2T+ed_LT_12+ed_gt_12+multiparous+
             smoke_house_outside, data=SalivaPart)
summary(modtest)

moddhea<-lm(meanlogdhea~male+BMI_2T+ed_LT_12+ed_gt_12+multiparous+
              smoke_house_outside, data=SalivaPart)
summary(moddhea)  #HOUSEHOLD SMOKE AND MALE SEX ASSOCIATED W/ DHEA

moddhea2<-lm(meanlogdhea~male+BMI_2T+ed_LT_12+ed_gt_12+multiparous+
              exsmoker, data=SalivaPart)
summary(moddhea2)

modmel<-lm(meanlogmel~male+BMI_2T+ed_LT_12+ed_gt_12+multiparous+
              smoke_house_outside, data=SalivaPart)
summary(modmel) # Borderline ass w/ melatonin and male sex

modcortslope<-lm(cortslope~male+BMI_2T+ed_LT_12+ed_gt_12+multiparous+
             smoke_house_outside+edad, data=SalivaPart)
summary(modcortslope) # Household smoke exposure assoicated with a larger cortisol slope

modcort1<-lm(newlogcort1~male+BMI_2T+ed_LT_12+ed_gt_12+multiparous+
                   smoke_house_outside+edad, data=SalivaPart)
summary(modcort1) #subtle 0.1 w/ value 1

modcort2<-lm(newlogcort2~male+BMI_2T+ed_LT_12+ed_gt_12+multiparous+
               smoke_house_outside, data=SalivaPart)
summary(modcort2) #Not associated with value 2 ).8

modcort3<-lm(newlogcort3~male+BMI_2T+ed_LT_12+ed_gt_12+multiparous+
               smoke_house_outside+edad, data=SalivaPart)
summary(modcort3) #0.1 at time 3 with a negative association

modProg2<-lm(meanlogprog~BPb__g_dL_M2T+male+BMI_2T+ed_LT_12+ed_gt_12+multiparous+
               smoke_house_outside, data=SalivaPart)
summary(modProg2)

modDHEA2<-lm(meanlogdhea~BPb__g_dL_M2T+male+BMI_2T+ed_LT_12+ed_gt_12+multiparous+
               smoke_house_outside, data=SalivaPart)
summary(modDHEA2)  #added lead with no difference

# NEED TO LOOK AT OTHER Hormones with lead
modest2<-lm(meanlogest~BPb__g_dL_M2T+male+BMI_2T+ed_LT_12+ed_gt_12+multiparous+
               smoke_house_outside, data=SalivaPart)
summary(modest2)

modtest2<-lm(meanlogtest~BPb__g_dL_M2T+male+BMI_2T+ed_LT_12+ed_gt_12+multiparous+
              smoke_house_outside, data=SalivaPart)
summary(modtest2)

modmel2<-lm(meanlogmel~BPb__g_dL_M2T+male+BMI_2T+ed_LT_12+ed_gt_12+multiparous+
               smoke_house_outside, data=SalivaPart)
summary(modmel2)

modcort2<-lm(cortslope~BPb__g_dL_M2T+male+BMI_2T+ed_LT_12+ed_gt_12+multiparous+
              smoke_house_outside, data=SalivaPart)
summary(modcort2)

# GESTATIONAL AGE WITH SPLINE

modspga<-gam(gest_age_days_d~s(meanlogprog,fx=TRUE,k=3) + multiparous+as.factor(ED_level)+s(BMI_2T, fx=TRUE,k=3)+
            smoke_house_outside, na.action=na.omit, data=SalivaPart)
summary(modspga)

plot(modspga,pch=18, scale=0, col=8, main ='natural spline', ylab="GestAge")

modspga2<-gam(gest_age_days_d~s(meanlogprog,fx=TRUE,k=3) + s(ahrr_mean, fx=TRUE, k=5)+ multiparous+as.factor(ED_level)+s(BMI_2T, fx=TRUE,k=3)+
               smoke_house_outside, na.action=na.omit, data=SP2)
summary(modspga2)

plot(modspga2,pch=18, scale=0, col=8, main ='natural spline', ylab="GestAge")

names(SalivaPart)

SP2$Cervix_PTGER2_Mean

