



#This is the LONG dataset

saliva$TubeNumber<-saliva$Tube.Number
saliva$TubeNumber

saliva$early1<-saliva$TubeNumber==1
saliva$early2<-saliva$TubeNumber==6
saliva$mid1<-saliva$TubeNumber==2
saliva$mid2<-saliva$TubeNumber==7
saliva$late1<-saliva$TubeNumber==5
saliva$late2<-saliva$TubeNumber==10
saliva$late3<-saliva$TubeNumber==4

saliva$time<-c(1, 2, 3)
saliva$time
#1, 2, 3 lined up properly with collection times
plot(saliva$time, saliva$early2)
plot(saliva$time, saliva$early1)
plot(saliva$time, saliva$mid1)
plot(saliva$time, saliva$mid2)
plot(saliva$time, saliva$late1)
plot(saliva$time, saliva$late2)
plot(saliva$time, saliva$late3)

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

saliva$logprog<-log(saliva$prog)
saliva$logcort<-log(saliva$cort)
saliva$logtest<-log(saliva$test)
saliva$logdhea<-log(saliva$dhea)
saliva$logest<-log(saliva$est)
saliva$logmel<-log(saliva$mel)

plot(saliva$logprog)
plot(salivaWide1$logprog1[saliva$logprog>0,])

salivaWide1<-saliva[saliva$time==1,]
salivaWide1$logprog1<-salivaWide1$logprog
salivaWide2<-saliva[saliva$time==2,]
salivaWide2$logprog2<-salivaWide2$logprog
salivaWide3<-saliva[saliva$time==3,]
salivaWide3$logprog3<-salivaWide3$logprog
plot(salivaWide$logprog2,salivaWide2$logprog )
plot(salivaWide$logprog3,salivaWide3$logprog )
plot(salivaWide$logprog1,salivaWide1$logprog )

salivaWide1$prog1<-salivaWide1$prog

salivaWide2$prog2<-salivaWide2$prog

salivaWide3$prog3<-salivaWide3$prog
plot(salivaWide3$prog3, salivaWide3$prog)
plot(salivaWide2$prog2, salivaWide2$prog)
plot(salivaWide1$prog1, salivaWide1$prog)


salivaWide<-cbind(salivaWide1, salivaWide2, salivaWide3)
names(salivaWide)


SW1<-salivaWide[, c("Folio", "prog1", "prog2", "prog3", "logprog1", "logprog2", "logprog3")]
#No Outliers- may need to bring some back in- goes to n=214 from 220
NMSW1<-na.omit(SW1)
NMSW2<-NMSW1[NMSW1$logprog1>0,]
NMSW3<-NMSW2[NMSW2$logprog2>0,]
NMSW4<-NMSW3[NMSW3$logprog3>0,]


cor(NMSW1$prog1, NMSW1$prog2)
cor(NMSW1$prog3, NMSW1$prog2)
cor(NMSW4$logprog1, NMSW4$logprog2)
cor(NMSW4$logprog3, NMSW4$logprog2)

#NEED TO IMPORT DATASET AND THEN Make box plots and then do bivariate

