AHRR_ALL <- read.csv("C:/R_Work_Dir/KforRv20141215.csv")

head(AHRR_ALL, 2 )
getwd()
KdfALL<-dput(names(AHRR_ALL))
KdfALL


AHRRDFALL <- AHRR_ALL[, c(mymethysites, mycovars)]
names(AHRRDFALL)
AHRRALLlong<-melt(AHRRDFALL, id.vars=c(mycovars)) #n=517 Similar findings, less significant for BWT/GA
head(AHRRALLlong)
#sensitivity analaysis if we had kept outliers
mod5ALL<-gee(value ~ BMI_2T + variable + gest_age_weeks_d + Fenton_Z_score +edad + 
             male+ ed_LT_12 +ed_gt_12 + smoke_house_outside+ Christiani_Lab+
             multiparous, id = folio, corstr = "exchangeable", data = AHRRALLlong)
summary(mod5ALL)$coefficients
print.model(mod5ALL)

ALLLM<-lm(ahrr_mean~BMI_2T  + gest_age_weeks_d + Fenton_Z_score +edad + 
     male+ ed_LT_12 +ed_gt_12 + smoke_house_outside+ Christiani_Lab+
     multiparous, data= AHRR_ALL)
summary(ALLLM)
d<-resid(ALLLM)
summary(d)
