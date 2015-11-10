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

# import the dataset
dat <- read_excel("H:/Epigenetics Bob Wright/Saliva/BPAAllan.xlsx", 1)
dat <- data.table(dat)
class(dat)

dat
# drop blank columns and rows
dat <- dat[!is.na(ID), names(dat)[-19], with = F]
names(dat)

# convert wide to long WAS ALREADY LONG
datmt <- melt(dat, id.vars = c("Run", "PDate", "Ptime", "Bdate", "Btime", "ID", 
                               "MCPP",  "MEP",	"MECPP",	"MEHHP",	"MBP",	"MiBP",	"MEOHP",	"MCMHP",	
                               "MBzP",	"MEHP",	"BPA"
))
datmt

# extract identifiers
datmt[, folio := as.numeric(gsub("MM-|-.T-M", "", ID))]
datmt[, etapa := gsub("MM-.{4}-|-M", "", ID)]
# what are the highest values
datmt[, max(value), by = variable]
datmt[which.max(value)]
dat[, describe(MBP)]

pthal2T<-datmt[etapa=="2T",]
pthal3T<-datmt[etapa=="3T",]

sal<-Saliva20150916

salphtal2Tc<-merge(sal,pthal2T,by="folio")
names(salphtal2Tc)

salphtal3Tc<-merge(sal,pthal3T,by="folio")
names(salphtal3Tc)


#"MCPP"                "MEP"                 "MECPP"               "MEHHP"              
#"MBP"                 "MiBP"                "MEOHP"               "MCMHP"              
#"MBzP"    "MEHP"                "BPA"   




#"MECPP"
#T2
fivenum(salphtal2Tc$MECPP)
summary(salphtal2Tc$MECPP)
sd(salphtal2Tc$MECPP)
summary(log(salphtal2Tc$MECPP))
sd(log(salphtal2Tc$MECPP))
t.test(log(MECPP)~preterm,data=salphtal2Tc)
exp(mean(log(salphtal2Tc$MECPP)))

MCPPlm<-lm(gest_age_days_d~log(MECPP), data=salphtal2Tc)
summary(MECPPlm)
resid(MECPPlm) #204
plot(MECPPlm)

mod_MECPP<-gam(gest_age_days_d~s(log(MECPP),fx=TRUE,k=3),
             na.action=na.omit, data=salphtal2Tc)
summary(mod_MECPP)
plot(mod_MECPP, pch=18, scale=0, col=12, main ='Second Trimster MECPP', ylab="Difference in Gestational Age (days)")

#T3
fivenum(salphtal3Tc$MECPP)
summary(salphtal3Tc$MECPP)
sd(salphtal3Tc$MECPP)
summary(log(salphtal3Tc$MECPP))
sd(log(salphtal3Tc$MECPP))
exp(mean(log(salphtal3Tc$MECPP)))
t.test(log(MECPP)~preterm,data=salphtal3Tc)

MECPPlm3<-lm(gest_age_days_d~log(MECPP), data=salphtal3Tc)
summary(MECPPlm3)
resid(MECPPlm3) #196
plot(MECPPlm3)

mod_MECPP3<-gam(gest_age_days_d~s(log(MECPP),fx=TRUE,k=3),
              na.action=na.omit, data=salphtal3Tc)
summary(mod_MECPP3)
plot(mod_MECPP3,pch=18, scale=0, col=12, main ='Third Trimster MECPP', ylab="Difference in Gestational Age (days)")


#MBP

#T2
fivenum(salphtal2Tc$MBP)
summary(salphtal2Tc$MBP)
sd(salphtal2Tc$MBP)
summary(log(salphtal2Tc$MBP))
sd(log(salphtal2Tc$MBP))
t.test(log(MBP)~preterm,data=salphtal2Tc)
exp(mean(log(salphtal2Tc$MBP)))

MBPlm<-lm(gest_age_days_d~log(MBP), data=salphtal2Tc)
summary(MBPlm)
resid(MBPlm) #204
plot(MBPlm)

mod_MBP<-gam(gest_age_days_d~s(log(MBP),fx=TRUE,k=3),
             na.action=na.omit, data=salphtal2Tc)
summary(mod_MBP)
plot(mod_MBP,ch=18, scale=0, col=12, main ='Second Trimster MBP', ylab="Difference in Gestational Age (days)")

#T3
fivenum(salphtal3Tc$MBP)
summary(salphtal3Tc$MBP)
sd(salphtal3Tc$MBP)
summary(log(salphtal3Tc$MBP))
sd(log(salphtal3Tc$MBP))
exp(mean(log(salphtal3Tc$MBP)))
t.test(log(MBP)~preterm,data=salphtal3Tc)

MBPlm3<-lm(gest_age_days_d~log(MBP), data=salphtal3Tc)
summary(MBPlm3)
resid(MBPlm3) #196
plot(MBPlm3)

mod_MBP3<-gam(gest_age_days_d~s(log(MBP),fx=TRUE,k=3),
             na.action=na.omit, data=salphtal3Tc)
summary(mod_MBP3)
plot(mod_MBP3,ch=18, scale=0, col=12, main ='Third Trimster MBP', ylab="Difference in Gestational Age (days)")


plot(log(salphtal2Tc$MBP),log(salphtal2Tc$MCPP))

#MCPP-> Issue with infinity or missingness?
fivenum(salphtal3Tc$MCPP)
summary(salphtal3Tc$MCPP)
sd(noInf2$MCPP)
summary(log(noInf2$MCPP))
sd(log(noInf2$MCPP))
t.test(log(MCPP)~preterm,data=noInf2)

noInf2<-salphtal2Tc[log(salphtal2Tc$MCPP)>-4,]
summary(log(salphtal2Tc$MCPP))
MCPPlm<-lm(gest_age_days_d~log(MCPP), data=noInf2)
plot(log(salphtal2Tc$MCPP),salphtal2Tc$gest_age_days_d )
summary(MCPPlm)
t<-resid(MCPPlm) #204
summary(t)
t
plot(MCPPlm)
names(salphtal2Tc)

mod_MCPP<-gam(gest_age_days_d~s(log(MCPP),fx=TRUE,k=3),
             na.action=na.omit, data=noInf2)
summary(mod_MCPP)
plot(mod_MCPP,ch=18, scale=0, col=12, main ='Second Trimster MCPP', ylab="Difference in Gestational Age (days)")

noInf2$preterm
table(noInf2$preterm) # 20 preterm 180 term
#T3

noInf3<-salphtal3Tc[log(salphtal3Tc$MCPP)>-10,]
fivenum(noInf3$MCPP)
summary(noInf3$MCPP)
sd(noInf3$MCPP)
summary(log(noInf3$MCPP))
sd(log(noInf3$MCPP))
t.test(log(MCPP)~preterm,data=noInf3)
table(noInf3$preterm)
exp(mean(log(noInf3$MCPP)))


MCPPlm3<-lm(gest_age_days_d~log(MCPP), data=salphtal3Tc[log(salphtal3Tc$MCPP)>-4,])
summary(MCPPlm3)
resid(MCPPlm3) #196
plot(MCPPlm3)

mod_MCPP3<-gam(gest_age_days_d~s(log(MCPP),fx=TRUE,k=3),
              na.action=na.omit, data=noInf3)
summary(mod_MCPP3)
plot(mod_MCPP3,ch=18, scale=0, col=12, main ="Third Trimster MCPP", ylab="Difference in Gestational Age (days)")
?lm

#MBzP
#T2
fivenum(salphtal2Tc$MBzP)
summary(salphtal2Tc$MBzP)
sd(salphtal2Tc$MBzP)
summary(log(salphtal2Tc$MBzP))
sd(log(salphtal2Tc$MBzP))
t.test(log(MBzP)~preterm,data=salphtal2Tc)
exp(mean(log(salphtal2Tc$MBzP)))

MBzPlm<-lm(gest_age_days_d~log(MBzP), data=salphtal2Tc)
summary(MBzPlm)
resid(MBzPlm) #204
plot(MBzPlm)

mod_MBzP<-gam(gest_age_days_d~s(log(MBzP),fx=TRUE,k=3),
             na.action=na.omit, data=salphtal2Tc)
summary(mod_MBzP)
plot(mod_MBzP,ch=18, scale=0, col=12, main ="Second Trimster MBzP", ylab="Difference in Gestational Age (days)")

#T3
fivenum(salphtal3Tc$MBzP)
summary(salphtal3Tc$MBzP)
sd(salphtal3Tc$MBzP)
summary(log(salphtal3Tc$MBzP))
sd(log(salphtal3Tc$MBzP))
t.test(log(MBzP)~preterm,data=salphtal3Tc)
exp(mean(log(salphtal3Tc$MBzP)))
MBzPlm3<-lm(gest_age_days_d~log(MBzP), data=salphtal3Tc)
summary(MBzPlm3)
resid(MBzPlm3) #196
plot(MBzPlm3)

mod_MBzP3<-gam(gest_age_days_d~s(log(MBzP),fx=TRUE,k=3),
              na.action=na.omit, data=salphtal3Tc)
summary(mod_MBzP3)
plot(mod_MBzP3,ch=18, scale=0, col=12, main ="Third Trimster MBzP", ylab="Difference in Gestational Age (days)")

plot(log(salphtal2Tc$MBzP),log(salphtal2Tc$MBzP))

#MiBP
#T2
fivenum(salphtal2Tc$MiBP)
summary(salphtal2Tc$MiBP)
sd(salphtal2Tc$MiBP)
summary(log(salphtal2Tc$MiBP))
sd(log(salphtal2Tc$MiBP))
t.test(log(MiBP)~preterm,data=salphtal2Tc)
exp(mean(log(salphtal2Tc$MiBP)))

MiBPlm<-lm(gest_age_days_d~log(MiBP), data=salphtal2Tc)
summary(MiBPlm)
resid(MiBPlm) #204
plot(MiBPlm)

mod_MiBP<-gam(gest_age_days_d~s(log(MBzP),fx=TRUE,k=3),
              na.action=na.omit, data=salphtal2Tc)
summary(mod_MiBP)
plot(mod_MiBP,scale=0, col=12, main ='Second Trimster MiBP', ylab="Difference in Gestational Age (days)")

#T3
fivenum(salphtal3Tc$MiBP)
summary(salphtal3Tc$MiBP)
sd(salphtal3Tc$MiBP)
summary(log(salphtal3Tc$MiBP))
sd(log(salphtal3Tc$MiBP))
t.test(log(MiBP)~preterm,data=salphtal3Tc)
exp(mean(log(salphtal3Tc$MiBP)))
MiBPlm3<-lm(gest_age_days_d~log(MiBP), data=salphtal3Tc)
summary(MiBPlm3)
resid(MiBPlm3) #196
plot(MiBPlm3)

mod_MiBP3<-gam(gest_age_days_d~s(log(MiBP),fx=TRUE,k=3),
               na.action=na.omit, data=salphtal3Tc)
summary(mod_MiBP3)
plot(mod_MiBP3,scale=0, col=12, main ='Third Trimster MiBP', ylab="Difference in Gestational Age (days)")

plot(log(salphtal2Tc$MiBP),log(salphtal2Tc$MiBP))

#MEOHP
#T2
fivenum(salphtal2Tc$MEOHP)
summary(salphtal2Tc$MEOHP)
sd(salphtal2Tc$MEOHP)
summary(log(salphtal2Tc$MEOHP))
sd(log(salphtal2Tc$MEOHP))
t.test(log(MEOHP)~preterm,data=salphtal2Tc)
exp(mean(log(salphtal2Tc$MEOHP)))

MEOHPlm<-lm(gest_age_days_d~log(MEOHP), data=salphtal2Tc)
summary(MEOHPlm)
resid(MEOHPlm) #204
plot(MEOHPlm)

mod_MEOHP<-gam(gest_age_days_d~s(log(MEOHP),fx=TRUE,k=3),
              na.action=na.omit, data=salphtal2Tc)
summary(mod_MEOHP)
plot(mod_MEOHP, scale=0, col=12, main ='Second Trimster MEOHP', ylab="Difference in Gestational Age (days)")

#T3

fivenum(salphtal3Tc$MEOHP)
summary(salphtal3Tc$MEOHP)
sd(salphtal3Tc$MEOHP)
summary(log(salphtal3Tc$MEOHP))
sd(log(salphtal3Tc$MEOHP))
t.test(log(MEOHP)~preterm,data=salphtal3Tc)
exp(mean(log(salphtal3Tc$MEOHP)))




MEOHPlm3<-lm(gest_age_days_d~log(MEOHP), data=salphtal3Tc)
summary(MEOHPlm3)
resid(MEOHPlm3) #196
plot(MEOHPlm3)

mod_MEOHP3<-gam(gest_age_days_d~s(log(MEOHP),fx=TRUE,k=3),
               na.action=na.omit, data=salphtal3Tc)
summary(mod_MEOHP3)
plot(mod_MEOHP3,scale=0, col=12, main ='Third Trimster MEOHP', ylab="Difference in Gestational Age (days)")

plot(log(salphtal2Tc$MEOHP),log(salphtal2Tc$MEOHP))

#MCMHP
#T2
fivenum(salphtal2Tc$MCMHP)
summary(salphtal2Tc$MCMHP)
sd(salphtal2Tc$MCMHP)
summary(log(salphtal2Tc$MCMHP))
sd(log(salphtal2Tc$MCMHP))
t.test(log(MCMHP)~preterm,data=salphtal2Tc)
exp(mean(log(salphtal2Tc$MCMHP)))






MCMHPlm<-lm(gest_age_days_d~log(MCMHP), data=salphtal2Tc)
summary(MCMHPlm)
resid(MCMHPlm) #204
plot(MCMHPlm)

mod_MCMHP<-gam(gest_age_days_d~s(log(MCMHP),fx=TRUE,k=3),
               na.action=na.omit, data=salphtal2Tc)
summary(mod_MCMHP)
plot(mod_MCMHP,scale=0, col=12, main ='Second Trimster MCMHP', ylab="Difference in Gestational Age (days)")
#T3
fivenum(salphtal3Tc$MCMHP)
summary(salphtal3Tc$MCMHP)
sd(salphtal3Tc$MCMHP)
summary(log(salphtal3Tc$MCMHP))
sd(log(salphtal3Tc$MCMHP))
t.test(log(MCMHP)~preterm,data=salphtal3Tc)
exp(mean(log(salphtal3Tc$MCMHP)))


MCMHPlm3<-lm(gest_age_days_d~log(MCMHP), data=salphtal3Tc)
summary(MCMHPlm3)
resid(MCMHPlm3) #196
plot(MCMHPlm3)

mod_MCMHP3<-gam(gest_age_days_d~s(log(MCMHP),fx=TRUE,k=3),
                na.action=na.omit, data=salphtal3Tc)
summary(mod_MCMHP3)
plot(mod_MCMHP3, scale=0, col=12, main ='Third Trimster MCMHP', ylab="Difference in Gestational Age (days)")



#MEHHP
#T2
fivenum(salphtal2Tc$MEHHP)
summary(salphtal2Tc$MEHHP)
sd(salphtal2Tc$MEHHP)
summary(log(salphtal2Tc$MEHHP))
sd(log(salphtal2Tc$MEHHP))
exp(mean(log(salphtal2Tc$MEHHP)))

t.test(log(MEHHP)~preterm,data=salphtal2Tc)


MEHHPlm<-lm(gest_age_days_d~log(MEHHP), data=salphtal2Tc)
summary(MEHHPlm)
resid(MEHHPlm) #204
plot(MEHHPlm)

mod_MEHHP<-gam(gest_age_days_d~s(log(MEHHP),fx=TRUE,k=3),
               na.action=na.omit, data=salphtal2Tc)
summary(mod_MEHHP)
plot(mod_MEHHP,scale=0, col=12, main ='Second Trimster MEHHP', ylab="Difference in Gestational Age (days)")

#T3
fivenum(salphtal3Tc$MEHHP)
summary(salphtal3Tc$MEHHP)
sd(salphtal3Tc$MEHHP)
summary(log(salphtal3Tc$MEHHP))
sd(log(salphtal3Tc$MEHHP))
exp(mean(log(salphtal3Tc$MEHHP)))

t.test(log(MEHHP)~preterm,data=salphtal3Tc)

MEHHPlm3<-lm(gest_age_days_d~log(MEHHP), data=salphtal3Tc)
summary(MEHHPlm3)
resid(MEHHPlm3) #196
plot(MEHHPlm3)

mod_MEHHP3<-gam(gest_age_days_d~s(log(MEHHP),fx=TRUE,k=3),
                na.action=na.omit, data=salphtal3Tc)
summary(mod_MEHHP3)
plot(mod_MEHHP3,scale=0, col=12, main ='Third Trimster MEHHP', ylab="Difference in Gestational Age (days)")

#MEHP
#T2
fivenum(salphtal2Tc$MEHP)
summary(salphtal2Tc$MEHP)
sd(salphtal2Tc$MEHP)
summary(log(salphtal2Tc$MEHP))
sd(log(salphtal2Tc$MEHP))
exp(mean(log(salphtal2Tc$MEHP)))
t.test(log(MEHP)~preterm,data=salphtal2Tc)

MEHPlm<-lm(gest_age_days_d~log(MEHP), data=salphtal2Tc)
summary(MEHPlm)
resid(MEHPlm) #204
plot(MEHPlm)

mod_MEHP<-gam(gest_age_days_d~s(log(MEHP),fx=TRUE,k=3),
               na.action=na.omit, data=salphtal2Tc)
summary(mod_MEHP)
plot(mod_MEHP, scale=0, col=12, main ='Second Trimster MEHP', ylab="Difference in Gestational Age (days)")

#T3

fivenum(salphtal3Tc$MEHP)
summary(salphtal3Tc$MEHP)
sd(salphtal3Tc$MEHP)
summary(log(salphtal3Tc$MEHP))
sd(log(salphtal3Tc$MEHP))
exp(mean(log(salphtal3Tc$MEHP)))
t.test(log(MEHP)~preterm,data=salphtal3Tc)

MEHPlm3<-lm(gest_age_days_d~log(MEHP), data=salphtal3Tc)
summary(MEHPlm3)
resid(MEHPlm3) #196
plot(MEHPlm3)

mod_MEHP3<-gam(gest_age_days_d~s(log(MEHP),fx=TRUE,k=3),
                na.action=na.omit, data=salphtal3Tc)
summary(mod_MEHP3)
plot(mod_MEHP3,scale=0, col=12, main ='Third Trimster MEHP', ylab="Difference in Gestational Age (days)")


#BPA
#T2
b2<-salphtal2Tc$BPA
fivenum(salphtal2Tc$BPA)
summary(salphtal2Tc$BPA)
sd(salphtal2Tc$BPA)
summary(log(salphtal2Tc$BPA))
sd(log(salphtal2Tc$BPA))
exp(mean(log(salphtal2Tc$BPA)))

t.test(log(BPA)~preterm,data=salphtal2Tc)
?t.test

BPAlm<-lm(gest_age_days_d~log(BPA), data=salphtal2Tc)
summary(BPAlm)
resid(BPAlm) #204
plot(BPAlm)

mod_BPA<-gam(gest_age_days_d~s(log(BPA),fx=TRUE,k=3),
              na.action=na.omit, data=salphtal2Tc)
summary(mod_BPA)
plot(mod_BPA, scale=0, col=12, main ='Second Trimster BPA', ylab="Difference in Gestational Age (days)")

#T3
fivenum(salphtal3Tc$BPA)
summary(salphtal3Tc$BPA)
sd(salphtal3Tc$BPA)
summary(log(salphtal3Tc$BPA))
sd(log(salphtal3Tc$BPA))
exp(mean(log(salphtal3Tc$BPA)))
t.test(log(BPA)~preterm,data=salphtal3Tc)

BPAlm3<-lm(gest_age_days_d~log(BPA), data=salphtal3Tc)
summary(BPAlm3)
resid(BPAlm3) #196
plot(BPAlm3)

mod_BPA3<-gam(gest_age_days_d~s(log(BPA),fx=TRUE,k=3),
               na.action=na.omit, data=salphtal3Tc)
summary(mod_BPA3)

plot(mod_BPA3, pch=18, scale=0, col=12, main ='Third Trimster BPA', ylab="Difference in Gestational Age (days)")


#MEP
#T2
b2<-salphtal2Tc$MEP
fivenum(salphtal2Tc$MEP)
summary(salphtal2Tc$MEP)
sd(salphtal2Tc$MEP)
summary(log(salphtal2Tc$MEP))
sd(log(salphtal2Tc$MEP))
exp(mean(log(salphtal2Tc$MEP)))
t.test(log(MEP)~preterm,data=salphtal2Tc)
?t.test

MEPlm<-lm(gest_age_days_d~log(MEP), data=salphtal2Tc)
summary(MEPlm)
resid(MEPlm) #204
plot(MEPlm)

mod_MEP<-gam(gest_age_days_d~s(log(MEP),fx=TRUE,k=3),
             na.action=na.omit, data=salphtal2Tc)
summary(mod_MEP)
plot(mod_MEP,scale=0, col=12, main ='Second Trimster MEP', ylab="Difference in Gestational Age (days)")

#T3
fivenum(salphtal3Tc$MEP)
summary(salphtal3Tc$MEP)
sd(salphtal3Tc$MEP)
summary(log(salphtal3Tc$MEP))
sd(log(salphtal3Tc$MEP))
exp(mean(log(salphtal3Tc$MEP)))
exp(4.88)
t.test(log(MEP)~preterm,data=salphtal3Tc)

MEPlm3<-lm(gest_age_days_d~log(MEP), data=salphtal3Tc)
summary(MEPlm3)
resid(MEPlm3) #196
plot(MEPlm3)

mod_MEP3<-gam(gest_age_days_d~s(log(MEP),fx=TRUE,k=3),
              na.action=na.omit, data=salphtal3Tc)
summary(mod_MEP3)

plot(mod_MEP3, pch=18, scale=0, col=12, main ='Third Trimster MEP', ylab="Difference in Gestational Age (days)")

##########

#MEP Most promising both in TTEST and SPline/LM

names(salphtal2Tc)
MEPlma<-lm(gest_age_days_d~log(MEP)+newBMI2T+edad+smoke_house_outside+ed_LT_12+ed_gt_12  , data=salphtal2Tc)
summary(MEPlma)
hist(salphtal2Tc$gest_age_days_d)

MEP3lma<-lm(gest_age_days_d~log(MEP)+newBMI2T+edad+smoke_house_outside+ed_LT_12+ed_gt_12  , data=salphtal3Tc)
summary(MEP3lma)

################

names(GASub)
