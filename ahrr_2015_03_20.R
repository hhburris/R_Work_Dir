library("colorspace", lib.loc="~/R/win-library/3.1")
library(reshape2)
library(gee)
library(ggplot2)
library(foreign)
library(mgcv)
library(stats)
library(rmarkdown)
library(knitr)
print(sessionInfo(), locale=FALSE)

getwd()
# import the csv file in this same folder - assuming that we opened this in RStudio
# such that the working directory is .../Heathers trials and tribulations
AHRR_GEE <- read.csv("AHRR_GEE.csv")
head(AHRR_GEE, 2)

Kdf<-dput(names(AHRR_GEE)) #restricted in SAS to non-outliers >20 for ahrr_pos1, >40 for ahrr_pos2, >20 for ahrr_pos3
Kdf

# here we specify a subset of the variables that we will need to go forward for use in this analysis
mymethysites <- c("ahrr_pos1", "ahrr_pos2", "ahrr_pos3")
# covariates and identifiers (one observation per person)

mycovars <- c("folio", "edad", "prepregnBMI", "Fenton_Z_score", "cbc_leuc", "cbc_neutroabs", "Christiani_Lab", "Plate_No", "Plate_Loc", "ahrr_mean", "preterm", "gest_age_days_d", "gest_age_weeks_d", "SGA", "AGA", "LGA", "ed_LT_12", "ed_12", "ed_gt_12", "ED_level", "multiparous", "smoke_house_outside", "male", "BMI_2T", "obese_2t", "overweight_2T", "normalweight_2T", "BMI_3Cat", "BWT_GA_3Cat", "folic_d", "cbc_linfo", "cbc_mono", "cbc_monoabs", "cbc_eosi", "cbc_baso", "cbc_hemat", "neut_perc_derived", "gran_perc_derived", "ahrr_pos1_sd", "ahrr_pos2_sd","ahrr_pos3_sd", "ahrr_mean_sd", "nr3_ob_sd")
# restricting to just these variables
#Still having trouble going to next line
ahrrdf2 <- AHRR_GEE[, c(mymethysites, mycovars)]

# create a long dataset where all of the covars will be repeated
# and each methylation observation will be one row
ahrrlong2 <- melt(ahrrdf2, id.vars = c(mycovars))
head(ahrrlong2, 3)

# we add "variable" as an indicator for methylation site (allow different intercepts)
# "value" is the methylation measure for a given folio-site combination
mod2 <- gee(value ~ BMI_2T + variable, id = folio, corstr = "exchangeable", data = ahrrlong2) 
summary(mod2)$coefficients
# highly significant!!! (estimate for BMI_2T is more than 2x the robust SE)
# see my old code for the summary function to pull out CI and p-value

# checking if LM is similar
modx<-lm(ahrr_pos1~BMI_2T, data=AHRR_GEE)
summary(modx)
mody<-lm(ahrr_pos2~BMI_2T, data=AHRR_GEE)
summary(mody)
modz<-lm(ahrr_pos3~BMI_2T, data=AHRR_GEE)
summary(modz)
#MEAN
modzz<-lm(ahrr_mean~BMI_2T, data=AHRR_GEE)
summary(modzz)
#GEE model VERY SIMILAR WITH MUCH BETTER P VALUE, but still same estimate
print.model(mod2)



table(ahrrdf2$SGA)

# for another day - we will link this up to a remote github repo
# git remote add origin https://github.com/hhburris/Heathers-K-anaylsis.git
# git push -u origin master

mod3 <- gee(value ~ BMI_2T + variable + Christiani_Lab, id = folio, corstr = "exchangeable", data = ahrrlong2) 
summary(mod2)$coefficients


mod4<- gee(value ~ BMI_2T + variable + Christiani_Lab + edad, id = folio, corstr = "exchangeable", data = ahrrlong2) 

mod5<-gee(value ~ BMI_2T + variable + Christiani_Lab + edad + gest_age_days_d, id = folio, corstr = "exchangeable", data = ahrrlong2)

mod5<-gee(value ~ BMI_2T + variable + Christiani_Lab + edad + gest_age_days_d + ed_LT_12 + ed_gt_12, id = folio, corstr = "exchangeable", data = ahrrlong2)
summary(mod5)$coefficients

print.model <- function(model, exp = F, digits = 2){
  se.coef <- if(class(model)[1] == "gee"){
    summary(model)$coefficients[, "Robust S.E."]} else {se.coef(model)}
  if(exp == T){
    or.df <- data.frame(OR = exp(coef(model)),
                        lower = exp(coef(model) -1.96 * se.coef),
                        upper = exp(coef(model) +1.96 * se.coef),
                        pvalue = round(1.96*(1-pnorm(abs(coef(model)/se.coef))), digits + 2))
  } else(
    or.df <- data.frame(estimate = (coef(model)),
                        lower = (coef(model) -1.96 * se.coef),
                        upper = (coef(model) +1.96 * se.coef),
                        pvalue = round(1.96*(1-pnorm(abs(coef(model)/se.coef))), digits + 2))
    
  )
  print(or.df, digits = digits)
}

# print it out!
print.model(mod5)

# FIGURE 1

ggplot(ahrrlong2, aes(variable, value, group = folio)) + 
  geom_point(position = position_jitter(height = 0, width = 0.1), alpha = 0.4) + 
  geom_line(alpha=0.2)+
  theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=17))+
  theme(axis.title.y = element_text(face="bold", size=17)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("DNA methylation (%5mC)")+
  xlab("AHRR CpG Site")+
  scale_x_discrete(breaks=c("ahrr_pos1", "ahrr_pos2", "ahrr_pos3"), labels=c("CpG 1", "CpG 2", "CpG 3"))

# Add this if we want title 
#ggtitle("AHRR DNA Methylation")+
#annotate("text", x =3.3,y = 80, label = "A", size=12, face="bold")
ggplot(ahrrlong2, aes(variable, value, group = folio)) + 
  geom_point(position = position_jitter(height = 0, width = 0.1), alpha = 0.4) + 
  geom_line(alpha=0.2)+
  theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=17))+
  theme(axis.title.y = element_text(face="bold", size=17)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("DNA methylation (%5mC)")+
  xlab("AHRR CpG Site")+
  annotate("text", x =1.3,y = 92, label = "r=0.83", size=7, face="bold")+
  annotate("text", x =2.7,y = 92, label = "r=0.61", size=7, face="bold")+
  scale_x_discrete(breaks=c("ahrr_pos1", "ahrr_pos2", "ahrr_pos3"), labels=c("CpG 1", "CpG 2", "CpG 3"))+
  theme(text = element_text(size=18, face="bold"))


ggplot(ahrrlong2, aes(BMI_2T, value, color = variable)) + 
  geom_point() + geom_smooth(method = "lm") + theme_bw(13) + 
  ylab("DNA methylation (%5mC)") + scale_color_hue("Site", labels=c("CpG 1", "CpG 2","CpG 3")) + 
  xlab("Maternal BMI")+
  ggtitle("Associations of Maternal BMI and Offspring AHRR DNA Methylation")

ggplot(ahrrlong2, aes(BMI_2T, value, color = variable)) + 
  geom_point() + geom_smooth(method = "lm") + theme_bw(13)+ theme(legend.title=element_text(size=16, face="bold"))+ 
  ylab("AHRR DNA methylation (%5mC)") + scale_color_hue("Site", labels=c("CpG 1", "CpG 2","CpG 3")) + 
  xlab("Maternal BMI")+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=17))+
  theme(axis.title.y = element_text(face="bold", size=17)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  theme(text = element_text(size=18, face="bold"))
#title   ggtitle("Associations of Maternal BMI and Offspring AHRR DNA Methylation")+


ggplot(ahrrlong2, aes(gest_age_weeks_d, value, color = variable)) + 
  geom_point() + geom_smooth(method = "lm") + theme_bw(13)+ theme(legend.title=element_text(size=16, face="bold"))+ 
  ylab("AHRR DNA methylation (%5mC)") + scale_color_hue("Site", labels=c("CpG 1", "CpG 2","CpG 3")) + 
  xlab("Gestational Age (weeks)")+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=17))+
  theme(axis.title.y = element_text(face="bold", size=17)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  theme(text = element_text(size=18, face="bold"))
#title:   ggtitle("Associations of Gestational Age and Offspring AHRR DNA Methylation")+


ggplot(ahrrlong2[ahrrlong2$gest_age_weeks_d >25,], aes(gest_age_weeks_d, value, color = variable)) + 
  geom_point() + geom_smooth(method = "lm") + theme_bw(13)+ theme(legend.title=element_text(size=16, face="bold"))+ 
  ylab("AHRR DNA methylation (%5mC)") + scale_color_hue("Site", labels=c("CpG 1", "CpG 2","CpG 3")) + 
  xlab("Gestational Age (weeks)")+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=17))+
  theme(axis.title.y = element_text(face="bold", size=17)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  theme(text = element_text(size=18, face="bold"))
#title: ggtitle("Associations of Gestational Age and AHRR DNA Methylation")+

ggplot(ahrrlong2, aes(Fenton_Z_score, value, color = variable)) + 
  geom_point() + geom_smooth(method = "lm") + theme_bw(13)+ theme(legend.title=element_text(size=16, face="bold"))+ 
  ylab("AHRR DNA methylation (%5mC)") + scale_color_hue("Site", labels=c("CpG 1", "CpG 2","CpG 3")) + 
  xlab("BWT for GA (Z score))")+
  theme(plot.title = element_text(size=17, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=17))+
  theme(axis.title.y = element_text(face="bold", size=17)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  theme(text = element_text(size=18, face="bold"))

#title : ggtitle("Associations of Birth Weight-For-Gestational Age and Offspring AHRR DNA Methylation")









# Box plot


# Maternal BMI Categories
a<-ggplot(AHRR_GEE, aes(factor(BMI_3Cat),ahrr_pos1))+
  geom_boxplot(size=0.2)+theme_bw()+
  theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=20))+
  theme(axis.title.y = element_text(face="bold", size=20)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("Mean AHRR DNA methylation (%5mC)")+
  xlab("Maternal BMI Categories")+
  scale_x_discrete(breaks=c("1", "2", "3"), labels=c("Normal Weight", "Overweight", "Obese"))+
  stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  theme(text = element_text(size=18, face="bold"))

names(AHRR_GEE)

# GA categories
b<-ggplot(AHRR_GEE, aes(factor(preterm),ahrr_pos1))+
  geom_boxplot(size=0.2)+theme_bw()+
  theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=20))+
  theme(axis.title.y = element_text(face="bold", size=20)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("AHRR (CpG 1) DNA methylation (%5mC)")+
  xlab("Gestational Age Categories")+
  scale_x_discrete(breaks=c("1", "0"), labels=c("Preterm", "Full term"))+
  stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  theme(text = element_text(size=18, face="bold"))

# Getting rid of missings
mycovars4 <- c("folio", "ahrr_pos1", "BWT_GA_3Cat")
newBoxDF<-AHRR_GEE[,c("folio", "ahrr_pos1", "BWT_GA_3Cat")]
nomiss<-na.omit(newBoxDF)


nomiss$BWT_GA_3Cat
######Plot for DNA methylation and BWT GA Categoreis
c<-ggplot(nomiss, aes(factor(BWT_GA_3Cat),ahrr_pos1))+
  geom_boxplot(size=0.4)+theme_bw()+
  theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=20))+
  theme(axis.title.y = element_text(face="bold", size=20)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("AHRR (CpG 1) DNA methylation (%5mC)")+
  xlab("Birth Weight-for-Gestational Age Categories")+
  scale_x_discrete(breaks=c("1", "2", "3"), labels=c("SGA", "AGA", "LGA"))+
  stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  theme(text = element_text(size=18, face="bold"))


summary(AHRR_GEE$BWT_GA_3Cat[is.na=TRUE])
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
multiplot(a, b, c, cols=3)



# Making it better with three across

d<-ggplot(AHRR_GEE, aes(factor(BMI_3Cat),ahrr_pos1))+
  geom_boxplot(size=0.2)+theme_bw()+
  theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=20))+
  theme(axis.title.y = element_text(face="bold", size=20)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("AHRR DNA methylation (%5mC)")+
  xlab("       ")+
  annotate("text", x =3.3,y = 80, label = "A", size=12, face="bold")+
  scale_x_discrete(breaks=c("1", "2", "3"), labels=c("Normal Weight", "Overweight", "Obese"))+
  stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  theme(text = element_text(size=18, face="bold"))



# GA categories
e<-ggplot(AHRR_GEE, aes(factor(preterm),ahrr_pos1))+
  geom_boxplot(size=0.2)+theme_bw()+
  theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=20))+
  theme(axis.title.y = element_text(face="bold", size=20)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("         ")+
  xlab("       ")+
  annotate("text", x =2.3,y = 80, label = "B", size=12, face="bold")+
  scale_x_discrete(breaks=c("1", "0"), labels=c("Preterm", "Full term"))+
  stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  theme(text = element_text(size=18, face="bold"))



######Plot for DNA methylation and BWT GA Categoreis
f<-ggplot(nomiss, aes(factor(BWT_GA_3Cat),ahrr_pos1))+
  geom_boxplot(size=0.4)+theme_bw()+
  theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=20))+
  theme(axis.title.y = element_text(face="bold", size=20)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("       ")+
  xlab("       ")+
  annotate("text", x =3.3,y = 80, label = "C", size=12, face="bold")+
  scale_x_discrete(breaks=c("1", "2", "3"), labels=c("SGA", "AGA", "LGA"))+
  stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  theme(text = element_text(size=16, face="bold"))

# Three panels for 
multiplot(d,e,f, cols=3)

l
# FIgure 3
l<-ggplot(ahrrlong2, aes(BMI_2T, value, color = variable)) + 
  geom_point() + geom_smooth(method = "lm") + theme_bw(13)+ theme(legend.title=element_text(size=16, face="bold"))+ 
  theme(legend.position=c(.8,.2))+
 scale_color_hue("Site", labels=c("CpG 1", "CpG 2","CpG 3")) + 
  xlab("Maternal BMI")+ ylab("AHRR DNA methylation (%5mC)")+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=17))+
  theme(axis.title.y = element_text(face="bold", size=17)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  theme(text = element_text(size=18, face="bold"))+
  annotate("text", x =44,y = 95, label = "A", size=12, face="bold")
l

#title   ggtitle("Associations of Maternal BMI and Offspring AHRR DNA Methylation")+


#ggplot(ahrrlong2, aes(gest_age_weeks_d, value, color = variable)) + 
# geom_point() + geom_smooth(method = "lm") + theme_bw(13)+ theme(legend.title=element_text(size=16, face="bold"))+ 
# ylab("AHRR DNA methylation (%5mC)") + scale_color_hue("Site", labels=c("CpG 1", "CpG 2","CpG 3")) + 
#xlab("Gestational Age (weeks)")+
#theme(plot.title = element_text(size=20, face="bold"))+
#theme(axis.title.x = element_text(face="bold", size=17))+
#theme(axis.title.y = element_text(face="bold", size=17)) +   
#theme(legend.title = element_text(size=16, face="bold"))+
#theme(text = element_text(size=18, face="bold"))
#title:   ggtitle("Associations of Gestational Age and Offspring AHRR DNA Methylation")+


m<-ggplot(ahrrlong2[ahrrlong2$gest_age_weeks_d >25,], aes(gest_age_weeks_d, value, color = variable)) + 
  geom_point() + geom_smooth(method = "lm") + theme_bw(13)+ theme(legend.title=element_text(size=16, face="bold"))+ 
   scale_color_hue("Site", labels=c("CpG 1", "CpG 2","CpG 3")) + 
  xlab("Gestational Age (weeks)")+ ylab("AHRR DNA methylation (%5mC)") +
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=17))+
  theme(axis.title.y = element_text(face="bold", size=17)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  theme(text = element_text(size=18, face="bold"))+
  annotate("text", x =44,y = 95, label = "B", size=12, face="bold")+
  theme(legend.position=c(.2,.2))
m

#title: ggtitle("Associations of Gestational Age and AHRR DNA Methylation")+

n<-ggplot(ahrrlong2, aes(Fenton_Z_score, value, color = variable)) + 
  geom_point() + geom_smooth(method = "lm") + theme_bw(13)+ theme(legend.title=element_text(size=16, face="bold"))+ 
  scale_color_hue("Site", labels=c("CpG 1", "CpG 2","CpG 3")) + 
  xlab("BWT for GA (Z score))")+ ylab("AHRR DNA methylation (%5mC)") +
  theme(plot.title = element_text(size=17, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=17))+
  theme(axis.title.y = element_text(face="bold", size=17)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  theme(text = element_text(size=18, face="bold"))+
  annotate("text", x =3,y = 95, label = "C", size=12, face="bold")+
  theme(legend.position=c(.8,.2))
n

multiplot(l,m,n, cols=3)

# Attempt at 3 across
# FIgure 3
l<-ggplot(ahrrlong2, aes(BMI_2T, value, color = variable)) + 
  geom_point() + geom_smooth(method = "lm") + theme_bw(13)+ theme(legend.title=element_text(size=16, face="bold"))+ 
  theme(legend.position="none")+
  scale_color_hue("Site", labels=c("CpG 1", "CpG 2","CpG 3")) + 
  xlab("Maternal BMI")+ ylab("   ")+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=17))+
  theme(axis.title.y = element_text(face="bold", size=17)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  theme(text = element_text(size=18, face="bold"))+
  annotate("text", x =44,y = 95, label = "A", size=6, face="bold")
l

#title   ggtitle("Associations of Maternal BMI and Offspring AHRR DNA Methylation")+


#ggplot(ahrrlong2, aes(gest_age_weeks_d, value, color = variable)) + 
# geom_point() + geom_smooth(method = "lm") + theme_bw(13)+ theme(legend.title=element_text(size=16, face="bold"))+ 
# ylab("AHRR DNA methylation (%5mC)") + scale_color_hue("Site", labels=c("CpG 1", "CpG 2","CpG 3")) + 
#xlab("Gestational Age (weeks)")+
#theme(plot.title = element_text(size=20, face="bold"))+
#theme(axis.title.x = element_text(face="bold", size=17))+
#theme(axis.title.y = element_text(face="bold", size=17)) +   
#theme(legend.title = element_text(size=16, face="bold"))+
#theme(text = element_text(size=18, face="bold"))
#title:   ggtitle("Associations of Gestational Age and Offspring AHRR DNA Methylation")+


m<-ggplot(ahrrlong2[ahrrlong2$gest_age_weeks_d >25,], aes(gest_age_weeks_d, value, color = variable)) + 
  geom_point() + geom_smooth(method = "lm") + theme_bw(13)+ theme(legend.title=element_text(size=16, face="bold"))+ 
  scale_color_hue("Site", labels=c("CpG 1", "CpG 2","CpG 3")) + 
  xlab("Gestational Age (weeks)")+ ylab("AHRR DNA methylation (%5mC)") +
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=17))+
  theme(axis.title.y = element_text(face="bold", size=17)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  theme(text = element_text(size=18, face="bold"))+
  annotate("text", x =43.5,y = 95, label = "B", size=6, face="bold")+
  theme(legend.position=c(.2,.2))
m

#title: ggtitle("Associations of Gestational Age and AHRR DNA Methylation")+

n<-ggplot(ahrrlong2, aes(Fenton_Z_score, value, color = variable)) + 
  geom_point() + geom_smooth(method = "lm") + theme_bw(13)+ theme(legend.title=element_text(size=16, face="bold"))+ 
  scale_color_hue("Site", labels=c("CpG 1", "CpG 2","CpG 3")) + 
  xlab("BWT for GA (Z score))")+ ylab(" ") +
  theme(plot.title = element_text(size=17, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=17))+
  theme(axis.title.y = element_text(face="bold", size=17)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  theme(text = element_text(size=18, face="bold"))+
  annotate("text", x =3,y = 95, label = "C", size=6, face="bold")+
  theme(legend.position="none")
n


multiplot(l,m,n, rows=3)

#####   Moving On 

# Stats for paper

cor(AHRR_GEE$ahrr_pos1, AHRR_GEE$ahrr_pos2, method="pearson") 
cor(AHRR_GEE$ahrr_pos2, AHRR_GEE$ahrr_pos3, method="pearson") 
cor(AHRR_GEE$ahrr_pos1, AHRR_GEE$ahrr_pos3, method="pearson") 

#Age- Not associated
modA <- gee(value ~ edad + variable, id = folio, corstr = "exchangeable", data = ahrrlong2) 

summary(modA)$coefficients
print.model(modA)

# Education
modB <- gee(value ~ ed_LT_12+ ed_gt_12 + variable, id = folio, corstr = "exchangeable", data = ahrrlong2) 
summary(modB)$coefficients
print.model(modB)

modB1 <- gee(value ~ ed_LT_12 + variable, id = folio, corstr = "exchangeable", data = ahrrlong2) 
summary(modB1)$coefficients
print.model(modB1)

#Male vs. female
modC<-gee(value ~ male + variable, id = folio, corstr = "exchangeable", data = ahrrlong2) 
summary(modC)$coefficients
print.model(modC)

#Parity
modD<-gee(value ~ multiparous + variable, id = folio, corstr = "exchangeable", data = ahrrlong2) 
summary(modD)$coefficients
print.model(modD)

#Household 
modE<-gee(value ~ smoke_house_outside + variable, id = folio, corstr = "exchangeable", data = ahrrlong2) 
summary(modE)$coefficients
print.model(modE)

#gestational age
modF<-gee(value ~ gest_age_weeks_d + variable, id = folio, corstr = "exchangeable", data = ahrrlong2) 
summary(modF)$coefficients
print.model(modF)

#BWT GA
modG<-gee(value ~ Fenton_Z_score + variable, id = folio, corstr = "exchangeable", data = ahrrlong2) 
summary(modG)$coefficients
print.model(modG)


mod2 <- gee(value ~ BMI_2T + variable, id = folio, corstr = "exchangeable", data = ahrrlong2) 
summary(mod2)$coefficients
print.model(mod2)

mod5<-gee(value ~ BMI_2T + variable + Christiani_Lab + edad + gest_age_days_d + ed_LT_12 + ed_gt_12, id = folio, corstr = "exchangeable", data = ahrrlong2)
summary(mod5)$coefficients

mod5a<-gee(value ~ BMI_2T + variable + Christiani_Lab +gest_age_weeks_d + Fenton_Z_score +edad + 
             male+ ed_LT_12 +ed_gt_12 + smoke_house_outside+
             multiparous, id = folio, corstr = "exchangeable", data = ahrrlong2)
summary(mod5a)$coefficients
print.model(mod5a)

# Christiani lab did not confound THUS the below is final model
mod5b<-gee(value ~ BMI_2T + variable + gest_age_weeks_d + Fenton_Z_score +edad + 
             male+ ed_LT_12 +ed_gt_12 + smoke_house_outside+
             multiparous, id = folio, corstr = "exchangeable", data = ahrrlong2)
summary(mod5b)$coefficients
print.model(mod5b)

# Categorical BMI
names(ahrrlong2)
mod5c<-gee(value ~ obese_2t+ overweight_2T + variable + gest_age_weeks_d + Fenton_Z_score +edad + 
             male+ ed_LT_12 +ed_gt_12 + smoke_house_outside+
             multiparous, id = folio, corstr = "exchangeable", data = ahrrlong2)
summary(mod5c)$coefficients
print.model(mod5c)

names(ahrrlong2)
mod5d<-gee(value ~ BMI_2T+ variable + preterm + Fenton_Z_score +edad + 
             male+ ed_LT_12 +ed_gt_12 + smoke_house_outside+
             multiparous, id = folio, corstr = "exchangeable", data = ahrrlong2)
summary(mod5d)$coefficients
print.model(mod5d)


mod5e<-gee(value ~ BMI_2T+ variable + gest_age_weeks_d + SGA+ LGA +edad + 
             male+ ed_LT_12 +ed_gt_12 + smoke_house_outside+
             multiparous, id = folio, corstr = "exchangeable", data = ahrrlong2)
summary(mod5e)$coefficients
print.model(mod5e)


names(ahrrlong2)

# sensitivity analyses

ahrrlong2$cbc_mono

mod5f<-gee(value ~ BMI_2T + variable + cbc_mono+ gran_perc_derived+gest_age_weeks_d +
             Fenton_Z_score +edad + 
             male+ ed_LT_12 +ed_gt_12 + smoke_house_outside+
             multiparous, id = folio, corstr = "exchangeable", data = ahrrlong2)
summary(mod5f)$coefficients
print.model(mod5f)

#nonmissing complete case with cbc


mycovars6 <- c("folio","BMI_2T", "gest_age_weeks_d", "Fenton_Z_score","edad",  "cbc_mono", 
               "gran_perc_derived",
               "male", "ed_LT_12", "ed_gt_12", "smoke_house_outside",
               "multiparous", "value","variable","ahrr_mean")
cbcDF<-ahrrlong2[,c(mycovars6)]
nomissCbc<-na.omit(cbcDF)

mod5g<-gee(value ~ BMI_2T + variable +gest_age_weeks_d +
             Fenton_Z_score +edad + 
             male+ ed_LT_12 +ed_gt_12 + smoke_house_outside+
             multiparous, id = folio, corstr = "exchangeable", data = nomissCbc)
summary(mod5g)$coefficients
print.model(mod5g)


mod5h<-gee(value ~ BMI_2T + variable + cbc_mono+ gran_perc_derived+gest_age_weeks_d +
             Fenton_Z_score +edad + 
             male+ ed_LT_12 +ed_gt_12 + smoke_house_outside+
             multiparous, id = folio, corstr = "exchangeable", data = nomissCbc)
summary(mod5h)$coefficients
print.model(mod5h)







#########  Let's see what happens to GIT
print(sessionInfo(), locale=FALSE)
 
# Figure out n's by running linear model and counting resids
resid(modzz)
#512 in unadjusted lm

modzz1<-lm(ahrr_mean~BMI_2T+ gest_age_weeks_d +
             Fenton_Z_score +edad + 
             male+ ed_LT_12 +ed_gt_12 + smoke_house_outside+
             multiparous, (data=AHRR_GEE))

summary(modzz1)
zz1<-resid(modzz1)
as.data.frame(zz1)

######## 507 in fully adjusted model

########CBC n=405 model

modzz2<-lm(ahrr_mean~BMI_2T+ cbc_mono+ gran_perc_derived+gest_age_weeks_d +
             Fenton_Z_score +edad + 
             male+ ed_LT_12 +ed_gt_12 + smoke_house_outside+
             multiparous, (data=AHRR_GEE))

summary(modzz2)
zz2<-resid(modzz2)
as.data.frame(zz2)


############## Folic acid 

modzz3a<-lm(ahrr_mean~BMI_2T+ folic_d+ gest_age_weeks_d +
             Fenton_Z_score +edad + 
             male+ ed_LT_12 +ed_gt_12 + smoke_house_outside+
             multiparous, (data=AHRR_GEE))

summary(modzz3a)
zz3a<-resid(modzz3a)
as.data.frame(zz3a)

modzz3<-gee(value ~ folic_d+ variable, 
            id = folio, corstr = "exchangeable", data = ahrrlong2)
summary(modzz3)$coefficients
print.model(modzz3)

# extract subest where folic_d not NA
mycovars5 <- c("folio","BMI_2T", "folic_d", "gest_age_weeks_d", "Fenton_Z_score","edad",  
                 "male", "ed_LT_12", "ed_gt_12", "smoke_house_outside",
                 "multiparous", "value","variable","ahrr_mean")
folicDF<-ahrrlong2[,c(mycovars5)]
nomissFolic<-na.omit(folicDF)

modzz4<-gee(value ~ folic_d+ variable, 
            id = folio, corstr = "exchangeable", data = ahrrlong2)
summary(modzz4)$coefficients
print.model(modzz4)


modzz5<-gee(value ~ BMI_2T + variable +gest_age_weeks_d +
             Fenton_Z_score +edad + 
             male+ ed_LT_12 +ed_gt_12 + smoke_house_outside+
             multiparous, id = folio, corstr = "exchangeable", data = nomissFolic)
summary(modzz5)$coefficients
print.model(modzz5)

#folic acid intake

modzz6<-gee(value ~ BMI_2T +folic_d+ variable +gest_age_weeks_d +
             Fenton_Z_score +edad + 
              male+ ed_LT_12 +ed_gt_12 + smoke_house_outside+
              multiparous, id = folio, corstr = "exchangeable", data = nomissFolic)
summary(modzz6)$coefficients
print.model(modzz6)


modzz7<-lm(ahrr_mean~BMI_2T +folic_d +gest_age_weeks_d +
             Fenton_Z_score +edad + 
             male+ ed_LT_12 +ed_gt_12 + smoke_house_outside+
             multiparous, (nomissFolic))
# n=65 for folic acid non-missings
summary(modzz7)
zz7<-resid(modzz7)
as.data.frame(zz7)

# Laborotory storage batch effects



mod5a<-gee(value ~ BMI_2T + variable + Christiani_Lab +gest_age_weeks_d + Fenton_Z_score +edad + 
             male+ ed_LT_12 +ed_gt_12 + smoke_house_outside+
             multiparous, id = folio, corstr = "exchangeable", data = ahrrlong2)
summary(mod5a)$coefficients
print.model(mod5a)

# Christiani lab did not confound THUS the below is final model
mod5b<-gee(value ~ BMI_2T + variable + gest_age_weeks_d + Fenton_Z_score +edad + 
             male+ ed_LT_12 +ed_gt_12 + smoke_house_outside+
             multiparous, id = folio, corstr = "exchangeable", data = ahrrlong2)
summary(mod5b)$coefficients
print.model(mod5b)

ChristLab<-ahrrlong2[ahrrlong2$Christiani_Lab=="1",]

mod5i<-gee(value ~ BMI_2T + variable + gest_age_weeks_d + Fenton_Z_score +edad + 
             male+ ed_LT_12 +ed_gt_12 + smoke_house_outside+
             multiparous, id = folio, corstr = "exchangeable", data = ChristLab)
summary(mod5i)$coefficients
print.model(mod5i)

Channing<-ahrrlong2[ahrrlong2$Christiani_Lab=="0",]
mod5j<-gee(value ~ BMI_2T + variable + gest_age_weeks_d + Fenton_Z_score +edad + 
             male+ ed_LT_12 +ed_gt_12 + smoke_house_outside+
             multiparous, id = folio, corstr = "exchangeable", data = Channing)
summary(mod5j)$coefficients
print.model(mod5j)

modzz8<-lm(ahrr_mean~BMI_2T+ gest_age_weeks_d +
             Fenton_Z_score +edad + 
             male+ ed_LT_12 +ed_gt_12 + smoke_house_outside+
             multiparous, (data=ChristLab))

summary(modzz8)
zz8<-resid(modzz8)
as.data.frame(zz8)
#Christiani lab 266

modzz9<-lm(ahrr_mean~BMI_2T+ gest_age_weeks_d +
             Fenton_Z_score +edad + 
             male+ ed_LT_12 +ed_gt_12 + smoke_house_outside+
             multiparous, (data=Channing))

summary(modzz9)
zz9<-resid(modzz9)
as.data.frame(zz9)
#Channing 241

# Self-reported BMI

#prepregnBMI

mod5PP<-gee(value ~ prepregnBMI + variable + gest_age_weeks_d + Fenton_Z_score +edad + 
             male+ ed_LT_12 +ed_gt_12 + smoke_house_outside+
             multiparous, id = folio, corstr = "exchangeable", data = ahrrlong2)
summary(mod5PP)$coefficients
print.model(mod5PP)

mean(ahrrlong2$prepregnBMI)
mean(ahrrlong2$BMI_2T)

###########After meeting with Andrea

# get rid of missing bwt ga
mcovars10<-c("folio", "Fenton_Z_score", "BMI_2T", "gest_age_weeks_d")
nomiss3<-AHRR_GEE[,c(mcovars10)]
nomissBWTGA<-na.omit(nomiss3)

cor(nomissBWTGA$BMI_2T, nomissBWTGA$Fenton_Z_score, method="pearson") 
cor(nomissBWTGA$BMI_2T, nomissBWTGA$gest_age_weeks_d, method="pearson")
cor(nomissBWTGA$gest_age_weeks_d, nomissBWTGA$Fenton_Z_score, method="pearson") 

plot(nomissBWTGA$gest_age_weeks_d, nomissBWTGA$Fenton_Z_score)
plot(nomissBWTGA$BMI_2T, nomissBWTGA$Fenton_Z_score)
plot(nomissBWTGA$BMI_2T, nomissBWTGA$gest_age_weeks_d)
