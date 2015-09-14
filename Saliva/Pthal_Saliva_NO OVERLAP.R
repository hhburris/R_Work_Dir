library(data.table)
library(StatMatch)
library(readxl)
library(Hmisc)
library(ggplot2)
library(reshape2)

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
"MCPP",	"MEP",	"MECPP",	"MEHHP",	"MBP",	"MiBP",	"MEOHP",	"MCMHP",	
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
sal$MBP2T<-pthal2T$MBP


#MBP
salphtal2Tb<-merge( pthal2T,sal, by="folio")
plot(salphtal2T$MBP)
SPH2T<-salphtal2T
SPH2Tb<-salphtal2Tb

SPH2Tb$meanlogdhea
names(SPH2Tb)
SPH2T$logMBP<-log(SPH2T$MBP)
plot(SPH2T$logMBP)
plot(SPH2T$logMBP,SPH2T$gest_age_days_d,)
MBPlm<-lm(gest_age_days_d~logMBP, data=SPH2T)
summary(MBPlm)
resid(MBPlm) #204 individuals w/ 2T MBP

#get rid of missings
plot
SPH2T$meanlogprog
proglm<-lm(meanlogprog)

##################NO OVERLAP BPA/SALIVA ####################
