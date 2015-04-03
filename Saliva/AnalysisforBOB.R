library("sas7bdat", lib.loc="~/R/win-library/3.1")
library("foreign", lib.loc="C:/Program Files/R/R-3.1.3/library")
library("grid", lib.loc="C:/Program Files/R/R-3.1.3/library")
library("stats", lib.loc="C:/Program Files/R/R-3.1.3/library")
library("nlme", lib.loc="C:/Program Files/R/R-3.1.3/library")
library("ggplot2", lib.loc="~/R/win-library/3.1")
library("mgcv", lib.loc="C:/Program Files/R/R-3.1.3/library")



merge3 <- read.csv("H:/Epigenetics Bob Wright/Saliva/merge3.csv") #after fix April 3, 2015
View(merge3)
names(merge3)


merge3$logprog1<-log(merge3$prog1)
merge3$logprog2<-log(merge3$prog2)
merge3$logprog3<-log(merge3$prog3)

merge3$logcort1<-log(merge3$cort1)
merge3$logcort2<-log(merge3$cort2)
merge3$logcort3<-log(merge3$cort3)

merge3$logtest1<-log(merge3$test1)
merge3$logtest2<-log(merge3$test2)
merge3$logtest3<-log(merge3$test3)

merge3$logdhea1<-log(merge3$dhea1)
merge3$logdhea2<-log(merge3$dhea2)
merge3$logdhea3<-log(merge3$dhea3)

merge3$logest1<-log(merge3$est1)
merge3$logest2<-log(merge3$est2)
merge3$logest3<-log(merge3$est3)

merge3$logmel1<-log(merge3$mel1)
merge3$logmel2<-log(merge3$mel2)
merge3$logmel3<-log(merge3$mel3)

noOL1<-merge3[merge3$logprog1>2,]
noOL2<-merge3[merge3$logprog2>2,]
noOL3<-merge3[merge3$logprog3>2,]

noOL4<-merge3[c(merge3$logprog1>2 & merge3$logprog2>2 & merge3$logprog3>2),]


# Progesterone



#box
a<-ggplot(merge3[merge3$logprog1>2,], aes(c("Tube 1"),logprog1))+
  geom_boxplot(size=0.2)+theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=20))+
  theme(axis.title.y = element_text(face="bold", size=20)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("Log Progresterone")+
  xlab("    ")+
  stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  theme(text = element_text(size=18, face="bold"))
a

b<-ggplot(merge3[merge3$logprog2>2,], aes(c("Tube 2"),logprog2))+
  geom_boxplot(size=0.2)+theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=20))+
  theme(axis.title.y = element_text(face="bold", size=20)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("           ")+
  xlab("    ")+
  stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  theme(text = element_text(size=18, face="bold"))
b

c<-ggplot(merge3[merge3$logprog3>2,], aes(c("Tube 3"),logprog3))+
  geom_boxplot(size=0.2)+theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=20))+
  theme(axis.title.y = element_text(face="bold", size=20)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("          ")+
  xlab("    ")+
  stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  theme(text = element_text(size=18, face="bold"))
c

multiplot(a, b,c, cols=3)


# to calculate means without outliers

merge4 <- read.csv("H:/Epigenetics Bob Wright/Saliva/merge5.csv") #after fix April 3, 2015 PM
View(merge4)

plot(merge4$newlogprog1, merge4$newlogprog2)
plot(merge4$newlogprog2, merge4$newlogprog3)

mod1<-lm(gest_age_days_d~meanlogprog, data=merge4)
summary(mod1)
plot(mod1)
plot(merge4$gestage_comb, merge4$meanlogprog)

names(merge4)

mean(merge4$meanlogprog[merge4$preterm==1], na.rm=TRUE)
mean(merge4$meanlogprog[merge4$preterm!=1],na.rm = TRUE)
t.test(merge4$meanlogprog~ merge4$preterm, na.rm=TRUE)
table(merge4$preterm)

merge4$meanlogprog

myvars<-c("Folio", "LOQ_BPAM3T","LOQ_BPAM2T","BPAM2T", "BPAM3T", "meanlogprog", "meanlogest",  
          "meanlogtest", "meanlogdhea","meanlogmel", 
          "newlogcort1", "newlogcort2", "newlogcort3",
          "meanlogcort",  "PAP", "PAP_3cat","PAP_cat",
          "LINEa_Pos1", "LINEa_Pos2", "LINEa_Pos3","LINEa_Pos4", 
          "LINEb_Pos1", "LINEb_Pos2",
          "LINEb_Pos3", "LINEb_Pos4", "LINEave_Pos1",  
          "LINEave_Pos2", "LINEave_Pos3",                    
          "LINEave_Pos4","Cervix_a_PTGER2_Pos_1", "Cervix_b_PTGER2_Pos_1","Cervix_PTGER2_Mean",
          "Fenton_Z_score" , "Fenton_Percentile", "gest_age_days_d", "gestage_comb",
          "ahrr_pos1" , "ahrr_pos2" ,"ahrr_pos3","ahrr_mean", "AXL_pos1", "AXL_pos2"   ,     "AXL_mean",
          "NR3_OB_Pos1"  ,   "TNF_Pos1"   ,    "TNF_Pos2"                        
          , "TNF_Pos3"  , "TNF_Pos4"  ,   "TNF_mean"                        
          , "IL8_Pos1"   ,   "IL6_Pos1"  ,  "IL6_Pos2"                        
          , "IL6_mean"  , "preterm" ,  "gest_age_days_d"                 
          , "gest_age_weeks_d"  , "SGA" , "AGA"                             
          , "LGA"  , "ed_LT_12"  ,"ed_12"                           
          ,"ed_gt_12" , "ED_level" ,  "multiparous", 
          "BMI_2T", "Christiani_Lab",
          "male", "smoke_house_outside"  )
table(merge4$male)
myvars
smalldf<-merge4[,myvars]

mod2<-lm(Fenton_Z_score~meanlogprog, data=smalldf)
summary(mod2)

mod3<-lm(Fenton_Z_score~meanlogcort, data=smalldf)
summary(mod3)

mod4<-lm(Fenton_Z_score~meanlogest, data=smalldf)
summary(mod4)

mod5<-lm(Fenton_Z_score~meanlogtest, data=smalldf)
summary(mod5)

mod6<-lm(Fenton_Z_score~meanlogmel, data=smalldf)
summary(mod6)

mod7<-lm(Fenton_Z_score~meanlogdhea, data=smalldf)
summary(mod7)

# saliva as outcome
plot(log(merge4$BPAM2T),merge4$meanlogprog)
mod8<-lm(meanlogprog~log(BPAM2T), data=smalldf) #need to deal with infinities
plot(log(merge4$BPAM3T), merge4$meanlogprog)
summary(mod8)
mod9<-lm(meanlogprog~BPAM3T, data=smalldf)
summary(mod9)

plot(gam8)
table(smalldf$LOQ_BPAM2T)
plot(smalldf$BPAM2T)
plot(smalldf$LOQ_BPAM2T)

#saliva and cervix

plot(smalldf$Cervix_PTGER2_Mean,smalldf$meanlogprog)
modcerv<-lm(Cervix_PTGER2_Mean~meanlogprog+PAP, data=smalldf) #N=54
summary(modcerv)
NN<-resid(modcerv)
summary(NN)

plot(smalldf$Cervix_PTGER2_Mean,smalldf$meanlogprog)
modcerv<-lm(Cervix_PTGER2_Mean~meanlogprog+PAP, data=smalldf) #N=54
summary(modcerv)
NN<-resid(modcerv)
summary(NN)

# Cortisone and NR3C1

plot(smalldf$meanlogcort,smalldf$NR3_OB_Pos1)
modNR3<-lm(NR3_OB_Pos1~meanlogcort+male+Fenton_Z_score, data=smalldf[smalldf$NR3_OB_Pos1>20,]) #N173- P 0.13
summary(modNR3)
RR<-resid(modNR3)
summary(RR)




# Taking a break to make the plots

d<-ggplot(merge4, aes(c("Tube 1"), newlogmel1))+
  geom_boxplot(size=0.2)+theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=20))+
  theme(axis.title.y = element_text(face="bold", size=20)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("Log Melatonin")+
  xlab("    ")+
  stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  theme(text = element_text(size=18, face="bold"))
d



e<-ggplot(merge4, aes(c("Tube 2"),newlogmel2))+
  geom_boxplot(size=0.2)+theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=20))+
  theme(axis.title.y = element_text(face="bold", size=20)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("          ")+
  xlab("    ")+
  stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  theme(text = element_text(size=18, face="bold"))
e

f<-ggplot(merge4, aes(c("Tube 3"),newlogmel3))+
  geom_boxplot(size=0.2)+theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=20))+
  theme(axis.title.y = element_text(face="bold", size=20)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("              ")+
  xlab("    ")+
  stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  theme(text = element_text(size=18, face="bold"))
f

multiplot(d,e,f, cols=3)


g<-ggplot(merge4, aes(c("Tube 1"), newlogest1))+
  geom_boxplot(size=0.2)+theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=20))+
  theme(axis.title.y = element_text(face="bold", size=20)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("Log Estradiol")+
  xlab("    ")+
  stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  theme(text = element_text(size=18, face="bold"))
g



h<-ggplot(merge4, aes(c("Tube 2"),newlogest2))+
  geom_boxplot(size=0.2)+theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=20))+
  theme(axis.title.y = element_text(face="bold", size=20)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("          ")+
  xlab("    ")+
  stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  theme(text = element_text(size=18, face="bold"))
h

i<-ggplot(merge4, aes(c("Tube 3"),newlogest3))+
  geom_boxplot(size=0.2)+theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=20))+
  theme(axis.title.y = element_text(face="bold", size=20)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("              ")+
  xlab("    ")+
  stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  theme(text = element_text(size=18, face="bold"))
i

multiplot(g, h, i, cols=3)




j<-ggplot(merge4, aes(c("Tube 1"), newlogtest1))+
  geom_boxplot(size=0.2)+theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=20))+
  theme(axis.title.y = element_text(face="bold", size=20)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("Log Testosterone")+
  xlab("    ")+
  stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  theme(text = element_text(size=18, face="bold"))
j



k<-ggplot(merge4, aes(c("Tube 2"),newlogtest2))+
  geom_boxplot(size=0.2)+theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=20))+
  theme(axis.title.y = element_text(face="bold", size=20)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("          ")+
  xlab("    ")+
  stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  theme(text = element_text(size=18, face="bold"))
k

l<-ggplot(merge4, aes(c("Tube 3"),newlogtest3))+
  geom_boxplot(size=0.2)+theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=20))+
  theme(axis.title.y = element_text(face="bold", size=20)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("              ")+
  xlab("    ")+
  stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  theme(text = element_text(size=18, face="bold"))
l

multiplot(j, k, l,  cols=3)


m<-ggplot(merge4, aes(c("Tube 1"), newlogdhea1))+
  geom_boxplot(size=0.2)+theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=20))+
  theme(axis.title.y = element_text(face="bold", size=20)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("Log DHEA")+
  xlab("    ")+
  stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  theme(text = element_text(size=18, face="bold"))
m



n<-ggplot(merge4, aes(c("Tube 2"),newlogdhea2))+
  geom_boxplot(size=0.2)+theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=20))+
  theme(axis.title.y = element_text(face="bold", size=20)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("          ")+
  xlab("    ")+
  stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  theme(text = element_text(size=18, face="bold"))
n

o<-ggplot(merge4, aes(c("Tube 3"),newlogdhea3))+
  geom_boxplot(size=0.2)+theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=20))+
  theme(axis.title.y = element_text(face="bold", size=20)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("              ")+
  xlab("    ")+
  stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  theme(text = element_text(size=18, face="bold"))
o

multiplot(m, n, o,  cols=3)


p<-ggplot(merge4, aes(c("Tube 1"), newlogcort1))+
  geom_boxplot(size=0.2)+theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=20))+
  theme(axis.title.y = element_text(face="bold", size=20)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("Log Cortisone")+
  xlab("    ")+
  stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  theme(text = element_text(size=18, face="bold"))
p



q<-ggplot(merge4, aes(c("Tube 2"),newlogcort2))+
  geom_boxplot(size=0.2)+theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=20))+
  theme(axis.title.y = element_text(face="bold", size=20)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("          ")+
  xlab("    ")+
  stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  theme(text = element_text(size=18, face="bold"))
q

r<-ggplot(merge4, aes(c("Tube 3"),newlogcort3))+
  geom_boxplot(size=0.2)+theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=20))+
  theme(axis.title.y = element_text(face="bold", size=20)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("              ")+
  xlab("    ")+
  stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  theme(text = element_text(size=18, face="bold"))
r

multiplot(p, q, r,  cols=3)


# SPAGHETTIN WITHOUT OUTLIERS (NONDETECTS)

A<-ggplot(saliva[saliva$logProg>2,], aes(time, logProg, group = Folio)) + 
  geom_point(position = position_jitter(height = 0, width = 0.1), alpha = 0.4) + 
  geom_line(alpha=0.2)+
  theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=17))+
  theme(axis.title.y = element_text(face="bold", size=17)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("log progesterone")+
  xlab("collection tube")+
  scale_x_discrete(limit=c("1", "2", "3"), 
                   labels=c("First", "Second", "Third"))+
  ggtitle("Progesterone")
A
B<-ggplot(saliva[log(saliva$test)>0.5,], aes(time, log(test), group = Folio)) + 
  geom_point(position = position_jitter(height = 0, width = 0.1), alpha = 0.4) + 
  geom_line(alpha=0.2)+
  theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=17))+
  theme(axis.title.y = element_text(face="bold", size=17)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("Log testosterone")+
  xlab("collection tube")+
  scale_x_discrete(limit=c("1", "2", "3"), 
                   labels=c("First", "Second", "Third"))+
  ggtitle("Testosterone")
B
C<-ggplot(saliva[log(saliva$dhea)>0.5,], aes(time, log(dhea), group = Folio)) + 
  geom_point(position = position_jitter(height = 0, width = 0.1), alpha = 0.4) + 
  geom_line(alpha=0.2)+
  theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=17))+
  theme(axis.title.y = element_text(face="bold", size=17)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("Log dhea")+
  xlab("collection tube")+
  scale_x_discrete(limit=c("1", "2", "3"), 
                   labels=c("First", "Second", "Third"))+
  ggtitle("DHEA")
C

D<-ggplot(saliva[log(saliva$est)>0.5,], aes(time, log(est), group = Folio)) + 
  geom_point(position = position_jitter(height = 0, width = 0.1), alpha = 0.4) + 
  geom_line(alpha=0.2)+
  theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=17))+
  theme(axis.title.y = element_text(face="bold", size=17)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("Log estrogen")+
  xlab("collection tube")+
  scale_x_discrete(limit=c("1", "2", "3"), 
                   labels=c("First", "Second", "Third"))+
  ggtitle("Estradiol")
D

E<-ggplot(saliva[log(saliva$mel)>-2,], aes(time, log(mel), group = Folio)) + 
  geom_point(position = position_jitter(height = 0, width = 0.1), alpha = 0.4) + 
  geom_line(alpha=0.2)+
  theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=17))+
  theme(axis.title.y = element_text(face="bold", size=17)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("Log melatonin")+
  xlab("collection tube")+
  scale_x_discrete(limit=c("1", "2", "3"), 
                   labels=c("First", "Second", "Third"))+
  ggtitle("Melatonin")
E
F<-ggplot(saliva, aes(time, log(cort), group = Folio)) + 
  geom_point(position = position_jitter(height = 0, width = 0.1), alpha = 0.4) + 
  geom_line(alpha=0.2)+
  theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=17))+
  theme(axis.title.y = element_text(face="bold", size=17)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("Log cortisones")+
  xlab("collection tube")+
  scale_x_discrete(limit=c("1", "2", "3"), 
                   labels=c("First", "Second", "Third"))+
  ggtitle("Cortisone")
F
multiplot(A, B, C, D, E, F, cols=3)


#############Note to self: Make correlations, put together a very brief slide show
