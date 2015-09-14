


sal$salGA<-salGA$GA_saliva_weeks
sal$salGA

mod_GASal<-lm(salGA~meanlogprog, data=sal)
mod_GASal
plot(mod_GASal)

plot(sal$prog1)
plot(sal$meanlogprog)

GASAL1<-ggplot(sal, aes(x=salGA, y= meanlogprog)) +
  geom_point(shape=1)+
  ylab("Log Progesterone") +    xlab("Gestational Age in weeks at collection")+
  ggtitle("Progesterone by week of gestation")
GASAL1

GASAL2<-ggplot(sal, aes(x=salGA, y= meanlogest)) +
  geom_point(shape=1)+
  ylab("Log estradiol") +    xlab("Gestational Age in weeks at collection")+
  ggtitle("Estradiol by week of gestation")
GASAL2

GASAL3<-ggplot(sal, aes(x=salGA, y= meanlogtest)) +
  geom_point(shape=1)+
  ylab("Log testosterone") +    xlab("Gestational Age in weeks at collection")+
  ggtitle("Testosterone by week of gestation")
GASAL3

GASAL4<-ggplot(sal, aes(x=salGA, y= meanlogdhea)) +
  geom_point(shape=1)+
  ylab("Log DHEA") +    xlab("Gestational Age in weeks at collection")+
  ggtitle("DHEA by week of gestation")
GASAL4







### Gestional age at birth lm

mod_GA1<-lm(gest_age_weeks_d~meanlogprog+salGA, data=sal)
summary(mod_GA1)
plot(mod_GA1)

mod_GA1sp<-gam(gest_age_days_d~s(meanlogprog,fx=TRUE,k=3) + multiparous+as.factor(ED_level)+
                 s(newBMI2T, fx=TRUE, k=3)+salGA+ s(meanlogcort, fx=TRUE, k=3)+ s(meanlogdhea, fx=TRUE,k=3)
+               smoke_house_outside, na.action=na.omit, data=sal)
summary(mod_GA1sp)

plot(mod_GA1sp,pch=18, scale=0, col=8, main ='natural spline', ylab="GestAge")
names(sal)

mod_GA2sp<-gam(gest_age_days_d~s(meanlogprog,fx=TRUE,k=3) + multiparous+as.factor(ED_level)+
                 s(newBMI2T, fx=TRUE, k=3)+salGA+
                 smoke_house_outside, na.action=na.omit, data=sal)
summary(mod_GA2sp)


mod_GAsAp<-gam(gest_age_days_d~s(meanlogprog,fx=TRUE,k=3) + 
                 s(meanlogdhea, fx=TRUE,k=3)+
                 s(meanlogcort, fx=TRUE,k=3)+
                 s(meanlogest, fx=TRUE,k=3)+
                 s(meanlogtest, fx=TRUE,k=3)+
                 s(meanlogmel, fx=TRUE,k=3)+
                 salGA,
                 na.action=na.omit, data=sal)

summary(mod_GAsAp)
plot(mod_GAsAp,pch=18, scale=0, col=12, main ='natural spline', ylab="GestAge")


mod_GADHEA<-gam(gest_age_days_d~
                 s(meanlogdhea, fx=TRUE,k=3)+
                 s(meanlogest, fx=TRUE,k=3)+
                 salGA+
                  edad+
                multiparous+as.factor(ED_level)+
                s(newBMI2T, fx=TRUE, k=3),
               na.action=na.omit, data=sal)

summary(mod_GADHEA)

plot(mod_GADHEA,pch=18, scale=0, col=10, main ='natural spline', ylab="GestAge")


mod_GADHEA2<-gam(gest_age_days_d~
                  meanlogdhea+
                  meanlogest+
                  salGA+
                  +as.factor(ED_level)+
                  s(newBMI2T, fx=TRUE, k=3),
                na.action=na.omit, data=sal)

summary(mod_GADHEA2)

plot(mod_GADHEA2,pch=18, scale=0, col=10, main ='natural spline', ylab="GestAge")

DHEA_ESTa<-ggplot(sal, aes(x=meanlogest, y=meanlogdhea)) +
  geom_point(shape=1,col=20, scale=2)+ ylab("Mean log-transformed Estradiol") +    xlab("Mean log-transformed DHEA")+
  ggtitle("Tube 1, DHEA and Estradiol, r=0.40")
DHEA_ESTa


# FOR PRESENTATION

mod_GAsAp<-gam(gest_age_days_d~s(meanlogprog,fx=TRUE,k=3) + 
                 s(meanlogdhea, fx=TRUE,k=3)+
                 s(meanlogcort, fx=TRUE,k=3)+
                 s(meanlogest, fx=TRUE,k=3)+
                 s(meanlogtest, fx=TRUE,k=3)+
                 s(meanlogmel, fx=TRUE,k=3)+
                 salGA,
               na.action=na.omit, data=sal)

summary(mod_GAsAp)
plot(mod_GAsAp,pch=18, scale=0, col=12, main ='natural spline', ylab="GestAge")

mod_PROG<-gam(gest_age_days_d~s(meanlogprog,fx=TRUE,k=3) + 
                 salGA,
               na.action=na.omit, data=sal)
summary(mod_PROG)
plot(mod_PROG,pch=18, scale=0, col=12, main ='Spline of hormone-gestational age relationship, adjusted for GA of collection', ylab="GestAge", xlab="Mean Log-transformed Progesterone")

mod_MEL<-gam(gest_age_days_d~s(meanlogmel,fx=TRUE,k=3) + 
                salGA,
              na.action=na.omit, data=sal)
summary(mod_MEL)
plot(mod_MEL,pch=18, scale=0, col=12, main ='Spline of hormone-gestational age relationship, adjusted for GA of collection', ylab="GestAge", xlab="Mean Log-transformed Melatonin")


mod_TEST<-gam(gest_age_days_d~s(meanlogtest,fx=TRUE,k=3) + 
               salGA,
             na.action=na.omit, data=sal)
summary(mod_TEST)
plot(mod_TEST,pch=18, scale=0, col=12, main ='Spline of hormone-gestational age relationship, adjusted for GA of collection',
     ylab="GestAge", xlab="Mean Log-transformed Testosterone")


mod_EST<-gam(gest_age_days_d~s(meanlogest,fx=TRUE,k=3) + 
                salGA,
              na.action=na.omit, data=sal)
summary(mod_EST)
plot(mod_EST,pch=18, scale=0, col=12, main ='Spline of hormone-gestational age relationship, adjusted for GA of collection',
     ylab="GestAge", xlab="Mean Log-transformed Estradiol")


mod_dhea<-gam(gest_age_days_d~s(meanlogdhea,fx=TRUE,k=3) + 
               salGA,
             na.action=na.omit, data=sal)
summary(mod_dhea)
plot(mod_dhea,pch=18, scale=0, col=12, main ='Spline of hormone-gestational age relationship, adjusted for GA of collection',
     ylab="GestAge", xlab="Mean Log-transformed DHEA")


# Now I'm adding EST to deal with dessication concerns (most highly correlated)

mod_dhea1<-gam(gest_age_days_d~
                s(meanlogdhea,fx=TRUE,k=3) + 
                s(meanlogest,fx=TRUE,k=3) +
                salGA,
              na.action=na.omit, data=sal)
summary(mod_dhea1)
plot(mod_dhea1,pch=18, scale=0, col=12, main ='Spline of hormone-gestational age relationship, adjusted for GA of collection')


#Adjust for multiple covariates
mod_dhea2<-gam(gest_age_days_d~
                 meanlogdhea + 
                 meanlogest +
                 s(newBMI2T, fx=TRUE, k=3)+
                 s(log(re_PbM2T), fx=TRUE, k=3)+
                 edad+
                 multiparous+
                 smoke_house_outside+
                 as.factor(ED_level)+
                 salGA,
               na.action=na.omit, data=sal)
summary(mod_dhea2)
plot(mod_dhea2,pch=18, scale=0, col=12, main ='Spline of BMI gestational age relationship, adjusted for GA of collection')

# More parsimonious
mod_dhea3<-gam(gest_age_days_d~
                 meanlogdhea + 
                 meanlogest +
                 smoke_house_outside+
                 ed_LT_12+
                 obese_2t+
                 overweight_2T+
                 salGA,
            
               na.action=na.omit, data=sal)
summary(mod_dhea3)
plot(mod_dhea3,pch=18, scale=0, col=12, main ='Spline of BMI gestational age relationship, adjusted for GA of collection')


# DHEA AS OUTCOME!
mod_DheaOUTa<-gam(meanlogdhea~s(log(re_PbM2T), fx=TRUE, k=3)+smoke_house_outside+edad+
                   ed_LT_12+salGA, na.action=na.omit, data=sal)
summary(mod_DheaOUTa) 
plot(mod_DheaOUTa,pch=18, scale=0, col=12, main ='Spline of DHEA-Pb relationship')

#to get ci

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



