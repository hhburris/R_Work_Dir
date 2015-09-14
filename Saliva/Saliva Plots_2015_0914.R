



# September 14, 2015

sal<-Saliva.for.R_2015_09_14

names(sal)


mysalvar<- c("folio", "meanlogprog","meanlogest", "meanlogtest", "meanlogdhea",
               "meanlogmel", "meanlogcort","newlogprog1"    ,     "newlogprog2"  ,       "newlogprog3" ,       
               "newlogcort1",   "newlogcort2"       ,  "newlogcort3"        , "newlogtest1"  ,      
              "newlogtest2"   ,      "newlogtest3" ,        "newlogest1"    ,      "newlogest2"   ,      
               "newlogest3"   ,       "newlogmel1"  ,        "newlogmel2"   ,       "newlogmel3",
                "newlogdhea1" ,        "newlogdhea2" ,        "newlogdhea3")

mysalvar

# Just the hormones and folio
sal2<-sal[,mysalvar]



# get rid of missings
nomisssal2<-na.omit(sal2)
plot(sal2$newlogprog1, sal2$newlogest1)
?ggplot
ggplot(nomisssal2, aes(x=newlogprog1, y=newlogprog2)) +
  geom_point(shape=2)
dev.off()

ggplot(sal2, aes(x=newlogprog1, y= newlogest1)) +
  geom_point(shape=1)+
ylab("Estradiol") +    xlab("Progesterone")+
  ggtitle("Tube 1, Estradiol and Progesterone, r=0.53")

ggplot(sal2, aes(x=newlogprog1, y= newlogcort1)) +
  geom_point(shape=1)+
  ylab("Cortisone") +    xlab("Progesterone")+
  ggtitle("Tube 1, Cortisone and Progesterone, r=0.60")

library(Hmisc)
rcorr(as.matrix(sal2[,mysalvarNoFolio]))

      mysalvarNoFolio<- c( "meanlogprog","meanlogest", "meanlogtest", "meanlogdhea",
                   "meanlogmel", "meanlogcort","newlogprog1"    ,     "newlogprog2"  ,       "newlogprog3" ,       
                   "newlogcort1",   "newlogcort2"       ,  "newlogcort3"        , "newlogtest1"  ,      
                   "newlogtest2"   ,      "newlogtest3" ,        "newlogest1"    ,      "newlogest2"   ,      
                   "newlogest3"   ,       "newlogmel1"  ,        "newlogmel2"   ,       "newlogmel3",
                   "newlogdhea1" ,        "newlogdhea2" ,        "newlogdhea3")
      



sal2<-sal[,mysalvar]

newBoxDFz<-merge4[,c("Folio", "meanlogprog", "preterm")]

nomissz<-na.omit(newBoxDFz)

plot(sal2)