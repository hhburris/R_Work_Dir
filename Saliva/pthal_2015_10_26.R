library(data.table)
library(StatMatch)
library(readxl)
library(Hmisc)
library(ggplot2)
library(reshape2)
library(mgcv)

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
#"MBzP"                "MEHP"                "BPA"        

#MBP

#T2
fivenum(salphtal2Tc$MBP)
summary(salphtal2Tc$MBP)
sd(salphtal2Tc$MBP)
summary(log(salphtal2Tc$MBP))
sd(log(salphtal2Tc$MBP))
t.test(log(MBP)~preterm,data=salphtal2Tc)

MBPlm<-lm(gest_age_days_d~log(MBP), data=salphtal2Tc)
summary(MBPlm)
resid(MBPlm) #204
plot(MBPlm)

mod_MBP<-gam(gest_age_days_d~s(log(MBP),fx=TRUE,k=3),
             na.action=na.omit, data=salphtal2Tc)
summary(mod_MBP)
plot(mod_MBP)

#T3
fivenum(salphtal3Tc$MBP)
summary(salphtal3Tc$MBP)
sd(salphtal3Tc$MBP)
summary(log(salphtal3Tc$MBP))
sd(log(salphtal3Tc$MBP))
t.test(log(MBP)~preterm,data=salphtal3Tc)

MBPlm3<-lm(gest_age_days_d~log(MBP), data=salphtal3Tc)
summary(MBPlm3)
resid(MBPlm3) #196
plot(MBPlm3)

mod_MBP3<-gam(gest_age_days_d~s(log(MBP),fx=TRUE,k=3),
             na.action=na.omit, data=salphtal3Tc)
summary(mod_MBP3)
plot(mod_MBP3)

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
plot(mod_MCPP)
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

MCPPlm3<-lm(gest_age_days_d~log(MCPP), data=salphtal3Tc[log(salphtal3Tc$MCPP)>-4,])
summary(MCPPlm3)
resid(MCPPlm3) #196
plot(MCPPlm3)

mod_MCPP3<-gam(gest_age_days_d~s(log(MCPP),fx=TRUE,k=3),
              na.action=na.omit, data=noInf3)
summary(mod_MCPP3)
plot(mod_MCPP3)
?lm

#MBzP
#T2
MBzPlm<-lm(gest_age_days_d~log(MBzP), data=salphtal2Tc)
summary(MBzPlm)
resid(MBzPlm) #204
plot(MBzPlm)

mod_MBzP<-gam(gest_age_days_d~s(log(MBzP),fx=TRUE,k=3),
             na.action=na.omit, data=salphtal2Tc)
summary(mod_MBzP)
plot(mod_MBzP)

#T3
MBzPlm3<-lm(gest_age_days_d~log(MBzP), data=salphtal3Tc)
summary(MBzPlm3)
resid(MBzPlm3) #196
plot(MBzPlm3)

mod_MBzP3<-gam(gest_age_days_d~s(log(MBzP),fx=TRUE,k=3),
              na.action=na.omit, data=salphtal3Tc)
summary(mod_MBzP3)
plot(mod_MBzP3)

plot(log(salphtal2Tc$MBzP),log(salphtal2Tc$MBzP))

#MiBP
#T2
MiBPlm<-lm(gest_age_days_d~log(MiBP), data=salphtal2Tc)
summary(MiBPlm)
resid(MiBPlm) #204
plot(MiBPlm)

mod_MiBP<-gam(gest_age_days_d~s(log(MBzP),fx=TRUE,k=3),
              na.action=na.omit, data=salphtal2Tc)
summary(mod_MiBP)
plot(mod_MiBP)

#T3
MiBPlm3<-lm(gest_age_days_d~log(MiBP), data=salphtal3Tc)
summary(MiBPlm3)
resid(MiBPlm3) #196
plot(MiBPlm3)

mod_MiBP3<-gam(gest_age_days_d~s(log(MiBP),fx=TRUE,k=3),
               na.action=na.omit, data=salphtal3Tc)
summary(mod_MiBP3)
plot(mod_MiBP3)

plot(log(salphtal2Tc$MiBP),log(salphtal2Tc$MiBP))

#MEOHP
#T2
MEOHPlm<-lm(gest_age_days_d~log(MEOHP), data=salphtal2Tc)
summary(MEOHPlm)
resid(MEOHPlm) #204
plot(MEOHPlm)

mod_MEOHP<-gam(gest_age_days_d~s(log(MEOHP),fx=TRUE,k=3),
              na.action=na.omit, data=salphtal2Tc)
summary(mod_MEOHP)
plot(mod_MEOHP)

#T3
MEOHPlm3<-lm(gest_age_days_d~log(MEOHP), data=salphtal3Tc)
summary(MEOHPlm3)
resid(MEOHPlm3) #196
plot(MEOHPlm3)

mod_MEOHP3<-gam(gest_age_days_d~s(log(MEOHP),fx=TRUE,k=3),
               na.action=na.omit, data=salphtal3Tc)
summary(mod_MEOHP3)
plot(mod_MEOHP3)

plot(log(salphtal2Tc$MEOHP),log(salphtal2Tc$MEOHP))

#MCMHP
#T2
MCMHPlm<-lm(gest_age_days_d~log(MCMHP), data=salphtal2Tc)
summary(MCMHPlm)
resid(MCMHPlm) #204
plot(MCMHPlm)

mod_MCMHP<-gam(gest_age_days_d~s(log(MCMHP),fx=TRUE,k=3),
               na.action=na.omit, data=salphtal2Tc)
summary(mod_MCMHP)
plot(mod_MCMHP)

#T3
MCMHPlm3<-lm(gest_age_days_d~log(MCMHP), data=salphtal3Tc)
summary(MCMHPlm3)
resid(MCMHPlm3) #196
plot(MCMHPlm3)

mod_MCMHP3<-gam(gest_age_days_d~s(log(MCMHP),fx=TRUE,k=3),
                na.action=na.omit, data=salphtal3Tc)
summary(mod_MCMHP3)
plot(mod_MCMHP3)



#MEHHP
#T2
MEHHPlm<-lm(gest_age_days_d~log(MEHHP), data=salphtal2Tc)
summary(MEHHPlm)
resid(MEHHPlm) #204
plot(MEHHPlm)

mod_MEHHP<-gam(gest_age_days_d~s(log(MEHHP),fx=TRUE,k=3),
               na.action=na.omit, data=salphtal2Tc)
summary(mod_MEHHP)
plot(mod_MEHHP)

#T3
MEHHPlm3<-lm(gest_age_days_d~log(MEHHP), data=salphtal3Tc)
summary(MEHHPlm3)
resid(MEHHPlm3) #196
plot(MEHHPlm3)

mod_MEHHP3<-gam(gest_age_days_d~s(log(MEHHP),fx=TRUE,k=3),
                na.action=na.omit, data=salphtal3Tc)
summary(mod_MEHHP3)
plot(mod_MEHHP3)

#MEHP
#T2
MEHPlm<-lm(gest_age_days_d~log(MEHP), data=salphtal2Tc)
summary(MEHPlm)
resid(MEHPlm) #204
plot(MEHPlm)

mod_MEHP<-gam(gest_age_days_d~s(log(MEHP),fx=TRUE,k=3),
               na.action=na.omit, data=salphtal2Tc)
summary(mod_MEHP)
plot(mod_MEHP)

#T3
MEHPlm3<-lm(gest_age_days_d~log(MEHP), data=salphtal3Tc)
summary(MEHPlm3)
resid(MEHPlm3) #196
plot(MEHPlm3)

mod_MEHP3<-gam(gest_age_days_d~s(log(MEHP),fx=TRUE,k=3),
                na.action=na.omit, data=salphtal3Tc)
summary(mod_MEHP3)
plot(mod_MEHP3)


#BPA
#T2
b2<-salphtal2Tc$BPA
fivenum(salphtal2Tc$BPA)
summary(salphtal2Tc$BPA)
sd(salphtal2Tc$BPA)
summary(log(salphtal2Tc$BPA))
sd(log(salphtal2Tc$BPA))

t.test(log(BPA)~preterm,data=salphtal2Tc)
?t.test

BPAlm<-lm(gest_age_days_d~log(BPA), data=salphtal2Tc)
summary(BPAlm)
resid(BPAlm) #204
plot(BPAlm)

mod_BPA<-gam(gest_age_days_d~s(log(BPA),fx=TRUE,k=3),
              na.action=na.omit, data=salphtal2Tc)
summary(mod_BPA)
plot(mod_BPA)

#T3
fivenum(salphtal3Tc$BPA)
summary(salphtal3Tc$BPA)
sd(salphtal3Tc$BPA)
summary(log(salphtal3Tc$BPA))
sd(log(salphtal3Tc$BPA))
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
plot(mod_MEP)

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
