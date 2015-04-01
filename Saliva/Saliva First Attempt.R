saliva<-Saliva_for_R
names(saliva)
saliva$cort<-saliva$Cortisone..ng.ml
saliva$cort
saliva$test<-saliva$Testosterone...pg.ml.
saliva$test
saliva$prog<-saliva$Progesterone...pg.ml.
saliva$prog
saliva$dhea<-saliva$DHEA...pg.ml.
saliva$dhea
saliva$est<-saliva$Estradiol..pg.ml
saliva$est
saliva$mel<-saliva$Melatonin...pg.ml.
saliva$mel

salivaWide1<-saliva[saliva$time==1,]
salivaWide1$LogProg1<-salivaWide1$logProg
salivaWide2<-saliva[saliva$time==2,]
salivaWide2$LogProg2<-salivaWide2$logProg
salivaWide3<-saliva[saliva$time==3,]
salivaWide3$LogProg3<-salivaWide3$logProg
plot(salivaWide$LogProg2,salivaWide2$logProg )
plot(salivaWide$LogProg3,salivaWide3$logProg )
plot(salivaWide$LogProg1,salivaWide1$logProg )

salivaWide1$prog1<-salivaWide1$prog

salivaWide2$prog2<-salivaWide2$prog

salivaWide3$prog3<-salivaWide3$prog
plot(salivaWide3$prog3, salivaWide3$prog)
plot(salivaWide2$prog2, salivaWide2$prog)
plot(salivaWide1$prog1, salivaWide1$prog)


salivaWide<-cbind(salivaWide1, salivaWide2, salivaWide3)
names(salivaWide)

SW1<-salivaWide[, c("Folio", "prog1", "prog2", "prog3", "LogProg1", "LogProg2", "LogProg3")]
#No Outliers- may need to bring some back in- goes to n=213 from 220
NMSW1<-na.omit(SW1)
NMSW2<-NMSW1[NMSW1$LogProg1>0,]
NMSW3<-NMSW2[NMSW2$LogProg2>0,]
NMSW4<-NMSW3[NMSW3$LogProg3>0,]


cor(NMSW1$prog1, NMSW1$prog2)
cor(NMSW1$prog3, NMSW1$prog2)
cor(NMSW4$LogProg1, NMSW4$LogProg2)
cor(NMSW4$LogProg3, NMSW4$LogProg2)

cov<-read.sas7bdat("cov_pyro.sas7bdat")
getwd()

SW2<-salivaWide
SW2$folio<-SW2$Folio
SW2$folio
total <- merge(cov,SW2,by="folio")
names(total)

summary(total$LogProg1)
hist(total$LogProg1)
hist(total$LogProg1[total$LogProg1>-Inf])
nomiss1<-total[total$LogProg1!=-Inf,]
nomiss2<-total[total$LogProg2!=-Inf,]
nomiss3<-total[total$LogProg3!=-Inf,]
mod1<-lm(gest_age_days_d~LogProg1, data=nomiss1)
summary(mod1)

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
print.model(mod1)
r<-resid(mod1)
r

mod2<-lm(gest_age_days_d~LogProg2, data=nomiss2)
summary(mod2)
noOL1<-nomiss1[nomiss1$LogProg1>0,]
gamfit<-gam(gest_age_days_d~s(LogProg1),data=noOL1)
summary(gamfit)
plot(gamfit)
plot(noOL1$LogProg1, noOL1$gest_age_days_d, main="Gest age",ylab="GestAgeDays")
abline
lines(noOL1$LogProg1, gamfit$fit,col="blue",lwd=3)

mod3<-lm(gest_age_days_d~LogProg3, data=nomiss3)
summary(mod3)

plot(total$LogProg1,total$gest_age_days_d )
plot(total$LogProg2,total$gest_age_days_d )
plot(total$LogProg3,total$gest_age_days_d )

# HIGHER IN ALL THREE AMONG PRETERM INFANTS 

mean(nomiss1$LogProg1[total$preterm==1])
mean(nomiss1$LogProg1[nomiss1$preterm!=1],na.rm = TRUE)
t.test(nomiss1$LogProg1~ nomiss1$preterm, na.rm=TRUE)
t.test(nomiss2$LogProg2~ nomiss2$preterm, na.rm=TRUE)
t.test(nomiss3$LogProg3~ nomiss3$preterm, na.rm=TRUE)

boxplot(nomiss1$LogProg1, nomiss1$preterm)
table(nomiss1$LogProg1, nomiss1$preterm)
table(nomiss1$preterm)

a<-ggplot(nomiss1, aes(factor(preterm),LogProg1))+
  geom_boxplot(size=0.2)+theme_bw()+
  theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=20))+
  theme(axis.title.y = element_text(face="bold", size=20)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("Log Progresterone 1")+
  xlab(" ")+
  scale_x_discrete(breaks=c("0", "1"), labels=c("Full term", "Preterm"))+
  stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  theme(text = element_text(size=18, face="bold"))
a

b<-ggplot(nomiss2, aes(factor(preterm),LogProg2))+
  geom_boxplot(size=0.2)+theme_bw()+
  theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=20))+
  theme(axis.title.y = element_text(face="bold", size=20)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("Log Progresterone 2")+
  xlab(" ")+
  scale_x_discrete(breaks=c("0", "1"), labels=c("Full term", "Preterm"))+
  stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  theme(text = element_text(size=18, face="bold"))
b
c<-ggplot(nomiss3, aes(factor(preterm),LogProg3))+
  geom_boxplot(size=0.2)+theme_bw()+
  theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=20))+
  theme(axis.title.y = element_text(face="bold", size=20)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  ylab("Log Progresterone 3")+
  xlab(" ")+
  scale_x_discrete(breaks=c("0", "1"), labels=c("Full term", "Preterm"))+
  stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  theme(text = element_text(size=18, face="bold"))
c

plot(salivaWide$LogProg3, salivaWide3$logProg)
salivaWide<-cbind(salivaWide1, salivaWide2, salivaWide3)
plot(salivaWide$LogProg1, salivaWide$LogProg2)
plot(noMissSalivaWide$LogProg2, noMissSalivaWide$LogProg3)
cor(noMissSalivaWide$LogProg1, noMissSalivaWide$LogProg2)
cor(noMissSalivaWide$LogProg2, noMissSalivaWide$LogProg3)
noMissSalivaWide<-na.omit(salivaWide)

?cbind
plot(salivaWide3$logProg, salivaWide2$logProg)
plot(salivaWide1$logProg, salivaWide2$logProg)
mod1<-lm(logProg~ logProg.2, na.omit, data=salivaWide)
names(salivaWide)

plot(saliva$cort)
plot(saliva$test)
plot(saliva$prog)
plot(saliva$dhea)
plot(saliva$est)
plot(saliva$mel)

saliva$TubeNumber<-saliva$Tube.Number
saliva$TubeNumber

saliva$early1<-saliva$TubeNumber==1
saliva$early2<-saliva$TubeNumber==6
saliva$mid1<-saliva$TubeNumber==2
saliva$mid2<-saliva$TubeNumber==7
saliva$late1<-saliva$TubeNumber==5
saliva$late2<-saliva$TubeNumber==10

saliva$time<-c(1, 2, 3)
saliva$time
#1, 2, 3 lined up properly with collection times
plot(saliva$time, saliva$early2)
plot(saliva$time, saliva$early1)
plot(saliva$time, saliva$mid1)
plot(saliva$time, saliva$mid2)
plot(saliva$time, saliva$late1)
plot(saliva$time, saliva$late2)


mean(saliva$prog,na.rm = TRUE)
fivenum(saliva$prog)
hist(saliva$prog)
saliva$logProg<-log(saliva$prog)
hist(saliva$logProg)

plot(saliva$logProg)
?corr
?cor
?anova
noOlLONG<-saliva[saliva$logProg>0,]
ggplot(noOlLONG, aes(Folio, logProg, color = as.factor(noOlLONG$time))) + 
  geom_point() + theme_bw(13)+ theme(legend.title=element_text(size=16, face="bold"))+ 
  ylab("Log Progesterone") + scale_color_hue("time", labels=c("first", "second","last")) + 
  xlab("Maternal BMI")+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.title.x = element_text(face="bold", size=17))+
  theme(axis.title.y = element_text(face="bold", size=17)) +   
  theme(legend.title = element_text(size=16, face="bold"))+
  theme(text = element_text(size=18, face="bold"))
#title   ggtitle("Associations of Maternal BMI and Offspring AHRR DNA Methylation")+
A1<-ggplot(noOlLONG, aes(time, logProg, group = Folio)) + 
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
ggtitle("Progesterone no outliers")
A1

A<-ggplot(saliva, aes(time, logProg, group = Folio)) + 
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

B<-ggplot(saliva, aes(time, log(test), group = Folio)) + 
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

C<-ggplot(saliva, aes(time, log(dhea), group = Folio)) + 
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

D<-ggplot(saliva, aes(time, log(est), group = Folio)) + 
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
ggtitle("Estrogen")

E<-ggplot(saliva, aes(time, log(mel), group = Folio)) + 
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
multiplot(A, B, C, D, E, F, cols=3 )
