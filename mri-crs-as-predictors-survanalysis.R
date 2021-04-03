library(survival)
library(survminer)


generate.Surv.obj<-function(allInclStudy){
ddsuba<-sub(" .*", "", corrected.data$Death.Date[allInclStudy])
odatea<-sub(" .*", "", corrected.data$OpDate[allInclStudy])
ddatea<-ifelse(ddsuba=="no" | ddsuba=="No", "2020-08-03", ddsuba)
# strptime(odate1, format = "%Y-%m-%d")
# strptime(ddate1, format = "%d/%m/%y")
os.timea<-difftime(ddatea, odatea, tz="GMT", units = "days")
#format="%dd%/%mm%/%yyyy%"
os.eventa<-as.numeric(ddatea!="2020-08-03")#is.deceased
#the relationship between any of the numeric variables and the log hazard is linear
return(Surv(os.timea,os.eventa))
}


check.Cox.assumpt<-function(){
allInclStudy<-intersect(allstudypop,which(!is.na(corrected.data$TumourVolume.mm.3.)))
length(allInclStudy)
allInclStudy2<-allInclStudy
agep<-corrected.data$AgeAtDiagnosis[allInclStudy]
cox.age<-coxph(Surv(os.timea,os.eventa)~agep)
plot(predict(cox.age),residuals(cox.age, type="martingale"),
     xlab="fitted values",ylab="Martingale residuals",
     main="Residual Plot")
plot(predict(cox.age),residuals(cox.age, type="martingale"),
     xlab="fitted values",ylab="Martingale residuals",
     main="Residual Plot",las=1)
#las=1 rotates the values on the Oy axis
#the values on the Oy axis are horizontallypositioned
#fitted=predicted values from cox.age numeric model
abline(h=0) # y=residual=0
lines(smooth.spline(predict(cox.age),residuals(cox.age,type="martingale")),col="red")

plot(predict(cox.age),residuals(cox.age, type="deviance"),
     xlab="fitted values",ylab="Deviance residuals",
     main="Residual Plot",las=1)
abline(h=0)
lines(smooth.spline(predict(cox.age),residuals(cox.age,type="deviance")),col="red")

#check linearity since the model contains numeric predictors
plot(predict(stratmodel3),residuals(stratmodel3, type="martingale"),
     xlab="fitted values",ylab="Martingale residuals",
     main="Residual Plot",las=1)
abline(h=0) # y=residual=0
lines(smooth.spline(predict(cox.age),residuals(cox.age,type="martingale")),col="red")

#cox.zph(coxph(Surv(os.timea,os.eventa)~stage.MRI.rescalef))

#will return a test for each of the individual variables or each of the coefficients, as well as
# a test for the overall model
cox.zph(cox.age)
#cox.mod<-coxph(Surv(os.timea,os.eventa)~stage.MRI.rescalef+agep)
#cox.zph(cox.mod)
#the test statistic is chisq
#GLOBAL=all variables at one
grade.histo<-factor(corrected.data$Grade[allInclStudy])
grade.histo<-grade.histo[-which(is.na(grade.histo))]
grade.histo.scale<-as.numeric(grade.histo)
grade.histo.rescale<-(grade.histo.scale-min(grade.histo.scale))/(max(grade.histo.scale)-min(grade.histo.scale))
sd(grade.histo.rescale)
summary(grade.histo.rescale)
grade.histo.rescalef<-cut(grade.histo.rescale,3)
levels(grade.histo.rescalef)<-levels(grade.histo)
compl.model<-coxph(Surv(os.timea,os.eventa)~stage.MRI+agep+tumour.vol)#+grade.histo)
cox.zph(compl.model)

#cox.zph(coxph(Surv(os.timea,os.eventa)~stage.MRI.rescalef+age.rescale+tumour.vol.rescale+factor(corrected.data$Grade[allInclStudy])))
early.st.pop<-intersect(allInclStudy,which(corrected.data$X.2 %in% c("IA","IB")))
non.adv.st.pop<-setdiff(allInclStudy,which(corrected.data$X.2 %in% c("IIIC1","IIIC2","IV")))
#if p>0.05=>we fail to reject the null hypothesis (that the hazards are proportional, 
#i.e. the hazard ratio is constant over time)
#we also build a plot to visualize this:
#if we allow the coefficients to change over time,
#or allow the hazard ratios to change over time,
#what changes we will see;
#if the coeff does not change over time,
#we would expect to see a change of 0
# a change of 0 means no change
ec.surv1<-generate.Surv.obj(early.st.pop)
stage.MRI1<-corrected.data$X.2[early.st.pop]
stage.MRI1<-factor(stage.MRI1)
age1<-corrected.data$AgeAtDiagnosis[early.st.pop]
tumour.vol1<-corrected.data$TumourVolume.mm.3.[early.st.pop]
grade.histo1<-corrected.data$Grade[early.st.pop]
grade.histo1<-factor(grade.histo1)
type.histo1<-corrected.data$X.1[early.st.pop]
type.histo1<-factor(type.histo1)
#chemo1<-chemotherapy[early.st.pop]
#radio1<-radiotherapy[early.st.pop]
coxph.model1<-coxph(ec.surv1~stage.MRI1+age1+tumour.vol1+grade.histo1+type.histo1)
par(mfrow=c(1,1))
plot(cox.zph(coxph.model1)[1])
abline(h=0,col="red")

par(mfrow=c(1,1))
plot(cox.zph(coxph.model1)[2])
abline(h=0,col="red")

par(mfrow=c(1,1))
plot(cox.zph(coxph.model1)[3])
abline(h=0,col="red")

par(mfrow=c(1,1))
plot(cox.zph(coxph.model1)[4])
abline(h=0,col="red")

summary(coxph.model1)
table(grade.histo1)
#only grade 2 sign; where is grade1
coxph.model2<-coxph(ec.surv1~stage.MRI1+age1+tumour.vol1+grade.histo1)
summary(coxph.model2)#concordance similar,but a bit less here
#4 explanatory variables, 319 observations, 57 events
#the limiting sample size is 57
#10*4=40; 15*4=60 (there are aprox. 15 observations per predictor variable)

#we can drop histological type without loss of predictive power
anova(coxph.model2,coxph.model1,test="LRT")#since models are nested
#pval large=>not a statistically significant difference between the two models
extractAIC(coxph.model1)
extractAIC(coxph.model2)#lower AIC for model2
#signoficantly better-fit model

coxph.model3<-coxph(ec.surv1~stage.MRI1+age1+tumour.vol1)
anova(coxph.model2,coxph.model3,test="LRT")
#p-val<<0.05=>significant difference between the models
extractAIC(coxph.model2)#lower AIC for model2
extractAIC(coxph.model3)
summary(coxph.model2)
summary(coxph.model3)#concordance is less

stage.histo1<-corrected.data$STAGE[early.st.pop]
stage.histo1<-factor(stage.histo1)
coxph.model4<-coxph(ec.surv1~stage.MRI1+age1+tumour.vol1+grade.histo1+stage.histo1)
anova(coxph.model2,coxph.model4,test="LRT")
#p-val<<0.05=>significant difference between the models
extractAIC(coxph.model2)
extractAIC(coxph.model4)#lower AIC for model4
summary(coxph.model2)
summary(coxph.model4)#larger, so better concordance
BIC(coxph.model2,coxph.model4)
AIC(coxph.model2,coxph.model4)
plot(survfit(coxph.model3,conf.int=FALSE),col="lightblue",xlab="OS Time (days)",ylab="Survival Probability",main="Comparison of Cox Regression Models")
lines(survfit(coxph.model2,conf.int=FALSE),col="red")
lines(survfit(coxph.model1,conf.int=FALSE),col="pink")
lines(survfit(coxph.model4,conf.int=FALSE),col="blue")
#library(aod)
#sigmap<-vcov(coxph.model4)
#wald.test(sigmap,Terms=1:5)
library(lmtest)
waldtest(coxph.model2,coxph.model4)
waldtest(coxph.model2)
waldtest(coxph.model1)
waldtest(coxph.model4)#even smaller pval
#The Wald test works by testing the null hypothesis that a set of parameters is equal to some value. 
#In the model being tested here, the null hypothesis is that 
#the 5 coefficients of interest are simultaneously equal to zero.
#pval very small=>reject the null hypothesis,
#this suggests that removing the variables from the model will 
#substantially harm the fit of that model
##it tests how far the estimated parameters are from zero 
#(or any other value under the null hypothesis) in standard errors
survdiff(ec.surv1~stage.MRI1)#+age1+tumour.vol1+grade.histo1+stage.histo1)
#The log rank test for difference in survival gives a p-value of p =
#we fail to reject the null hypothesis 
#these independent groups do not differ significantly in survival 

#we only have 5 independent variables that are tested simultaneously
#so use FDR?
surv.data<-data.frame(ec.surv1[,1],stage.MRI1,age1,tumour.vol1,grade.histo1,stage.histo1)
library(fuzzySim)
FDR(data=surv.data,sp.cols=1,var.cols=2:6)

library(jtools)
library(ggstance)
library(broom.mixed)
plot_summs(coxph.model3,coxph.model2,coxph.model1,coxph.model4)
#the farthest from zero, the more predictive 

#stage.MRI1+age1+tumour.vol1+grade.histo1+stage.histo1
summary(coxph(ec.surv1 ~ stage.MRI1))$coef #not significant
summary(coxph(ec.surv1 ~ age1))$coef
summary(coxph(ec.surv1 ~ tumour.vol1))$coef
summary(coxph(ec.surv1 ~ grade.histo1))$coef
summary(coxph(ec.surv1 ~ stage.histo1))$coef#not significant for grade2
summary(coxph(ec.surv1 ~ type.histo1))$coef#only sign for endometrioid

# library(survival)
# cox.age1<-coxph(ec.surv1~age1)
# plot(predict(cox.age1),residuals(cox.age1, type="martingale"),
#      xlab="fitted values",ylab="Martingale residuals",
#      main="Residual Plot")
# plot(predict(cox.age1),residuals(cox.age1, type="martingale"),
#      xlab="fitted values",ylab="Martingale residuals",
#      main="Residual Plot",las=1)
#las=1 rotates the values on the Oy axis
#the values on the Oy axis are horizontallypositioned
#fitted=predicted values from cox.age numeric model
# abline(h=0) # y=residual=0
# lines(smooth.spline(predict(cox.age1),residuals(cox.age1,type="martingale")),col="red")
# 
# plot(predict(cox.age1),residuals(cox.age1, type="deviance"),
#      xlab="fitted values",ylab="Deviance residuals",
#      main="Residual Plot",las=1)
# abline(h=0)
# lines(smooth.spline(predict(cox.age1),residuals(cox.age1,type="deviance")),col="red")


cox.agetum<-coxph(ec.surv1~age1+tumour.vol1)
# plot(predict(cox.agetum),residuals(cox.agetum, type="martingale"),
#      xlab="fitted values",ylab="Martingale residuals",
#      main="Residual Plot")
plot(predict(cox.agetum),residuals(cox.agetum, type="martingale"),
     xlab="fitted values",ylab="Martingale residuals",
     main="Residual Plot",las=1)
#las=1 rotates the values on the Oy axis
#the values on the Oy axis are horizontallypositioned
#fitted=predicted values from cox.age numeric model
abline(h=0) # y=residual=0
lines(smooth.spline(predict(cox.agetum),residuals(cox.agetum,type="martingale")),col="red")

plot(predict(cox.agetum),residuals(cox.agetum, type="deviance"),
     xlab="fitted values",ylab="Deviance residuals",
     main="Residual Plot",las=1)
abline(h=0)
lines(smooth.spline(predict(cox.agetum),residuals(cox.agetum,type="deviance")),col="red")

##for the whole model containing numerical variables
plot(predict(coxph.model4),residuals(coxph.model4, type="martingale"),
     xlab="fitted values",ylab="Martingale residuals",
     main="Residual Plot",las=1)
#las=1 rotates the values on the Oy axis
#the values on the Oy axis are horizontallypositioned
#fitted=predicted values from cox.age numeric model
abline(h=0) # y=residual=0
lines(smooth.spline(predict(coxph.model4),residuals(coxph.model4,type="martingale")),col="red")


early.st.dat<-corrected.data[early.st.pop,]
ggcoxfunctional(ec.surv1 ~ age1 + log(age1) + sqrt(age1),data=early.st.dat)


#the order of predictors#stage.MRI1+age1+tumour.vol1+grade.histo1+stage.histo1
coxph.model5<-coxph(ec.surv1~stage.MRI1+age1+tumour.vol1+strata(grade.histo1))#+strata(stage.histo1))

par(mfrow=c(1,1))
plot(cox.zph(coxph.model4)[1])
abline(h=0,col="red")

par(mfrow=c(1,1))
plot(cox.zph(coxph.model4)[2])#stratify for age
abline(h=0,col="red")

par(mfrow=c(1,1))
plot(cox.zph(coxph.model4)[3])#stratify for tumour volume
abline(h=0,col="red")

par(mfrow=c(1,1))
plot(cox.zph(coxph.model4)[4])#stratify for grade histo!
abline(h=0,col="red")

par(mfrow=c(1,1))
plot(cox.zph(coxph.model4)[5])#stratify for stage histo!
abline(h=0,col="red")


par(mfrow=c(1,1))
plot(cox.zph(coxph.model5)[1])
abline(h=0,col="red")

par(mfrow=c(1,1))
plot(cox.zph(coxph.model5)[2])#stratify for age
abline(h=0,col="red")

par(mfrow=c(1,1))
plot(cox.zph(coxph.model5)[3])#stratify for tumour volume
abline(h=0,col="red")
}


search4confounders<-function(){
plot(tumour.vol1,ec.surv1[,1])
plot(age1,ec.surv1[,1])
plot(survfit(ec.surv1 ~ stage.MRI1), col=c("black", "red"), fun="cloglog")
plot(survfit(ec.surv1 ~ grade.histo1), col=c("black", "red"), fun="cloglog")
plot(survfit(ec.surv1 ~ strata(grade.histo1)), col=c("black", "red"), fun="cloglog")
#stage.MRI1+age1+tumour.vol1+strata(grade.histo1)
}


stratify<-function(){
age.groups1<-vector()
age.groups1[which(age1<45)]<-'<45'
age.groups1[which(age1>=45&age1<=70)]<-'45-70'
age.groups1[which(age1>70)]<-'>70'
tumourv.groups1<-vector()
tumourv.groups1[which(tumour.vol1<100000)]<-'<0.0001m3'
tumourv.groups1[which(tumour.vol1>=100000)]<-'0.0001-0.0004m3'
tumourv.groups2<-vector()
tumourv.groups2[which(tumour.vol1<=median(tumour.vol1))]<-'<=8839.888mm3'
tumourv.groups2[which(tumour.vol1>median(tumour.vol1))]<-'>8839.888mm3'

nostratmodel<-coxph(ec.surv1~stage.MRI1+age1+tumour.vol1+grade.histo1)
stratmodel<-coxph(ec.surv1~stage.MRI1+age1+strata(tumourv.groups1)+grade.histo1)#+strata(stage.histo1))
ggsurvplot(survfit(nostratmodel), fun = "cloglog",data=surv.data[2:5])

#log-log plots
ggsurvplot(survfit(ec.surv1~tumourv.groups1),fun="cloglog",data=surv.data[2:5])
ggsurvplot(survfit(stratmodel),fun="cloglog", data = surv.data[2:5], conf.int = FALSE)

#KM-curves
ggsurvplot(survfit(stratmodel), data = surv.data[2:5], conf.int = FALSE)
ggsurvplot(survfit(ec.surv1~tumourv.groups1),data=surv.data[2:5])


stratmodel3<-coxph(ec.surv1~stage.MRI1+age1+strata(tumourv.groups2)+grade.histo1)
ggsurvplot(survfit(ec.surv1 ~ (tumour.vol1 > median(tumour.vol1, na.rm = T)), data = surv.data[2:5]), data = surv.data[2:5])
ggsurvplot(survfit(ec.surv1 ~ tumourv.groups2, data = surv.data[2:5]), data = surv.data[2:5])

#log-log plots
ggsurvplot(survfit(ec.surv1~tumourv.groups2),fun="cloglog",data=surv.data[2:5])
ggsurvplot(survfit(stratmodel3),fun="cloglog", data = surv.data[2:5], conf.int = FALSE)

#KM-curves
ggsurvplot(survfit(stratmodel3), data = surv.data[2:5], conf.int = FALSE)

cox.zph(nostratmodel)
cox.zph(stratmodel)
stratmodel2<-coxph(ec.surv1~stage.MRI1+age1+strata(tumourv.groups1,grade.histo1))
#stratmodel vs 2: compare concordance and p-vals (they are not nested)

cox.zph(stratmodel2)
cox.zph(stratmodel4)

stratmodel4<-coxph(ec.surv1~stage.MRI1+age1+strata(tumourv.groups2,grade.histo1))
#stratmodel3 vs 4: compare concordance and p-vals (they are not nested)

#stratmodel5<-coxph(ec.surv1~stage.MRI1+age.groups1+strata(tumourv.groups2,grade.histo1))
#model not good, it gives warning

#log-log plots
ggsurvplot(survfit(ec.surv1~grade.histo1),fun="cloglog",data=surv.data[2:5])
ggsurvplot(survfit(stratmodel2),fun="cloglog", data = surv.data[2:5], conf.int = FALSE)
ggsurvplot(survfit(stratmodel4),fun="cloglog", data = surv.data[2:5], conf.int = FALSE)

#KM-curves
ggsurvplot(survfit(stratmodel2), data = surv.data[2:5], conf.int = FALSE)
ggsurvplot(survfit(stratmodel4), data = surv.data[2:5], conf.int = FALSE)
ggsurvplot(survfit(ec.surv1~grade.histo1),data=surv.data[2:5])

#ggsurvplot(survfit(ec.surv1~stage.MRI1+age.groups1+strata(tumourv.groups2,grade.histo1)),data = surv.data[2:5], conf.int = FALSE)
ggsurvplot(survfit(stratmodel4), data=surv.data[2:5], #color = "#2E9FDF",
           ggtheme = theme_minimal())

par(mfrow=c(1,1))
plot(cox.zph(stratmodel)[1])
abline(h=0,col="red")
plot(cox.zph(stratmodel)[2])#stratify for age
abline(h=0,col="red")
plot(cox.zph(stratmodel)[3])#stratify for tumour volume
abline(h=0,col="red")

par(mfrow=c(1,1))
plot(cox.zph(stratmodel3)[1])
abline(h=0,col="red")
plot(cox.zph(stratmodel3)[2])#stratify for age
abline(h=0,col="red")
plot(cox.zph(stratmodel3)[3])#stratify for tumour volume
abline(h=0,col="red")

tumourv.groups1<-factor(tumourv.groups1)
tumourv.groups2<-factor(tumourv.groups2)
#levels(tumourv.groups1)<-c("<0.0001m3","0.0001-0.0004m3")
#log-log plots
splots1<-list()
#splots1[[1]]
g1<-ggsurvplot(survfit(ec.surv1~tumourv.groups1),xlab="Time on a logarithmic scale",legend.title = "tumour volume",legend.labs =c("< 100000 mm^3","\u2265 100000 mm^3"),fun="cloglog",data=surv.data[2:5],title="with tumour grouping I alone")
plot(survfit(ec.surv1~tumourv.groups1), fun = "cloglog", xlab = "Log(Time)",
     ylab = "Log-log survival") 
#splots1[[2]]
g2<-ggsurvplot(survfit(ec.surv1~tumourv.groups2),xlab="Time on a logarithmic scale", fun="cloglog",legend.title = "tumour volume",legend.labs = c("\u2264 8839.888 mm^3","> 8839.888 mm^3"),data=surv.data[2:5],title="with tumour grouping II alone")
plot(survfit(ec.surv1~tumourv.groups2), fun = "cloglog", xlab = "Log(Time)",
     ylab = "Log-log survival") 
#ggsurvplot(survfit(stratmodel),fun="cloglog", data = surv.data[2:5], conf.int = FALSE)
#ggsurvplot(survfit(stratmodel3),fun="cloglog", data = surv.data[2:5], conf.int = FALSE)
#labels(splots1[[1]],"A")
title<-grid.text("Complementary log-log plots for Survival Models",gp=gpar(fontsize=20,fontface="bold",fontfamily="Calibri"))
png("checkPHAssumption-tumour_volume.png",width=16,height=8,units='in',res=300)
ag<-grid.arrange(g1$plot,g2$plot,ncol=2)#arrange_ggsurvplots(splots1)#"Checking the PH assumption with complementary log-log plots\n")
botm<-grid.text("n = 319\n(subset: MRI stage = I)",gp=gpar(fontsize=16,fontfamily="Arial"))
grid.arrange(title,ag,layout_matrix=rbind(c(1,1),c(2,2),c(2,2),c(2,2),c(2,2)),bottom=botm)
dev.off()
}


compare.categorized.vars<-function(){
splots2=list()
#KM-curves
splots2[[1]]
splots21<-ggsurvplot(survfit(ec.surv1~tumourv.groups1),xlab="OS Time (days)",title="Kaplan-Meier estimate with tumour volume grouping I alone",legend.title = "tumour volume",legend.labs = c("< 100000 mm^3","\u2265 100000 mm^3"),data=surv.data[2:5],risk.table = TRUE)
splots2[[3]]
p21<-surv_pvalue(
  fit=survfit(ec.surv1~tumourv.groups1),
  data = surv.data[2:5],
  method = "survdiff",
  test.for.trend = FALSE,
  combine = TRUE
)
p21val<-p21$pval
splots21$plot<-splots21$plot+annotate("text", x = 500, y=0.2, label = paste("p = ",p21val,"\nfor log-rank test"),size=5)
splots22<-ggsurvplot(survfit(ec.surv1 ~ (tumour.vol1 > median(tumour.vol1, na.rm = T))),xlab="OS Time (days)",title="Kaplan-Meier estimate with tumour volume grouping II alone",legend.title = "tumour volume",legend.labs = c("\u2264 8839.888 mm^3","> 8839.888 mm^3"), data = surv.data[2:5],risk.table = TRUE)
#ggsurvplot(survfit(ec.surv1 ~ tumourv.groups2, data = surv.data[2:5]), data = surv.data[2:5])
splots2[[2]]
p22<-surv_pvalue(
  fit=survfit(ec.surv1 ~ (tumour.vol1 > median(tumour.vol1, na.rm = T))),
  data = surv.data[2:5],
  method = "survdiff",
  test.for.trend = FALSE,
  combine = TRUE
)
p22val<-p22$pval
splots22$plot<-splots22$plot+annotate("text", x = 1000, y=0.2, label = paste("p = ",p22val," for log-rank test"),size=5)
#stratmodel
splots23<-ggsurvplot(survfit(stratmodel),xlab="OS Time (days)",title="\nSurvival curves predicted by the Cox (PH) model with MRI stage adjusted\nfor age + strata(tumour.groupingI) + tumour grade",legend.title="tumour volume",legend.labs = c("< 100000 mm^3","\u2265 100000 mm^3"), data = surv.data[2:5], conf.int = FALSE,risk.table = TRUE)#, pval = TRUE)
fit23<-surv_fit(stratmodel$formula,data = surv.data[2:5])
p23<-surv_pvalue(
  fit=fit23,
  data = surv.data[2:5],
  method = "survdiff",
  test.for.trend = FALSE,
  combine = TRUE
)
p23val<-p23$pval
splots23$plot<-splots23$plot+annotate("text", x = 1000, y=0.2, label = paste("p = ",p23val," for log-rank test"),size=5)

splots2[[4]]
splots24<-ggsurvplot(survfit(stratmodel3),xlab="OS Time (days)",title="\nSurvival curves predicted by the Cox (PH) model with MRI stage adjusted\nfor age + strata(tumour.groupingII) + tumour grade", legend.labs = c("\u2264 8839.888 mm^3","> 8839.888 mm^3"),legend.title="tumour volume", data = surv.data[2:5], conf.int = FALSE,risk.table = TRUE)#, pval = TRUE)
fit24<-surv_fit(stratmodel3$formula,data = surv.data[2:5])
p24<-surv_pvalue(
  fit=fit24,
  data = surv.data[2:5],
  method = "survdiff",
  test.for.trend = FALSE,
  combine = TRUE
)
p24val<-p24$pval
splots24$plot<-splots24$plot+annotate("text", x = 1000, y=0.2, label = paste("p = ",p24val," for log-rank test"),size=5)
title<-grid.text("The Overall Survival Separated by Tumour Volume (in Grouping I vs II)",gp=gpar(fontsize=20,fontface="bold",fontfamily="Calibri"))#,fontfamily="Calibri"))
subtitle<-grid.text("Survival Curves for Model Comparisons\n",gp=gpar(fontsize=16,fontfamily="Arial"))
bottm<-grid.text("\nn = 319\n(subset: MRI stage = I)",gp=gpar(fontsize=16,fontfamily="Arial"))
png("compareTumourGroupingsWmodels.png",width=16,height=12,units='in',res=300)
#arrange_ggsurvplots(splots2,ncol=2,nrow=2,title=,labels=c("C","D","E","F"))
gra1<-grid.arrange(splots21$plot,splots22$plot,splots21$table,splots22$table,layout_matrix=rbind(c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(3,4),c(3,4)))
gra2<-grid.arrange(splots23$plot,splots24$plot,splots23$table,splots24$table,layout_matrix=rbind(c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(3,4),c(3,4)))
gra12<-grid.arrange(gra1,gra2,nrow=2)
grid.arrange(title,subtitle,gra12,
             heights = unit.c(grobHeight(title) + 1.2*margin, 
                              grobHeight(subtitle) + margin, 
                              #grobHeight(bot)+marginb,
                              unit(1,"null")), bottom=bottm)
dev.off()
#labels=c("C","D","E","F")
survdiff(ec.surv1~stage.MRI1+age1+strata(tumourv.groups2)+grade.histo1)
survdiff(ec.surv1~stage.MRI1+age1+strata(tumourv.groups1)+grade.histo1)
levels(tumourv.groups2)
}


model.wMRInCRS<-function(){
tumourvolstrat<-strata(tumourv.groups2)
stratifmodel<-coxph(ec.surv1~stage.MRI1+age1+tumourvolstrat+grade.histo1)
#plot
ggsurvplot(survfit(stratifmodel),data=surv.data[2:5], conf.int = TRUE, palette = "Dark2")#, 
           #censor = FALSE)
ggsurvplot(survfit(nostratmodel1),data=surv.df2, conf.int = TRUE, palette = "Dark2")
ggsurvplot(survfit(nostratmodel2),data=surv.df2, conf.int = TRUE, palette = "Dark2")
fit3 <- list(complete.with.both = survfit(nostratmodel1), without.MRI.stage = survfit(nostratmodelwout1),without.CRS = survfit(nostratmodelwout2), MRI.stage.better.fit=survfit(stratifmodel))
fin.df<-data.frame(ec.surv1[,1],ec.surv1[,2],stage.MRI1,age1,tumour.vol1,grade.histo1,clinicrs)
fin.df2<-data.frame(ec.surv1[,1],ec.surv1[,2],stage.MRI1,age1,tumourv.groups2,grade.histo1,clinicrs)
png("MRIstage&CRSinOS.png",width=16,height=8,units='in',res=300)
ggls<-list("with MRI stage + CRS + age + tumour.groupingII"=survfit(model11nostrat),"with CRS + age + tumour.groupingII"=survfit(model11a),"with MRI stage + age + tumour.groupingII"=survfit(model11bc))
ggsc<-ggsurvplot_combine(ggls,xlab="OS Time (days)",title="Comparing the Primary Predictors Strength (MRI stage vs CRS): Cox (PH) Regression Model Comparisons by Predicted Survival Curves",legend.title = "Model by predictors", fin.df2)#, palette = "Dark2")
ggsc$plot<-ggsc$plot+scale_colour_manual(values = c(rgb(1,0,0,1),rgb(0,0,1,0.5),rgb(0.9,0.65,0)))
#legend.labs =c("with MRI stage and CRS","without MRI stage","without CRS","with MRI stage and histo grade"),
#Comparison between estimated survival curves for different predictors of OS: MRI stage vs CRS
grid.arrange(ggsc$plot,bottom="\nn = 319, number of events = 57 (subset: MRI stage = I)")
dev.off()
}


plotID.predstrength.survCurvComp<-function(subpop){
  ggls<-list("with MRI stage + CRS + age + tumour.groupingII"=survfit(model11nostrat,subset=subpop),"with CRS + age + tumour.groupingII"=survfit(model11a,subset=subpop),"with MRI stage + age + tumour.groupingII"=survfit(model11bc,subset=subpop))
  ggsc<-ggsurvplot_combine(ggls,xlab="OS Time (days)",title="Comparing the Primary Predictors Strength (MRI stage vs CRS): Cox (PH) Regression Model Comparisons by Survival Curves Estimate",legend.title = "Model by predictors", fin.df2[subpop,])#, palette = "Dark2")
  ggsc$plot<-ggsc$plot+scale_colour_manual(values = c(rgb(1,0,0,1),rgb(0,0,1,0.5),rgb(0.9,0.65,0)))
  return(ggsc$plot)
}

plotM.predstrength.survCurvComp<-function(model11nostrat,model11a,model11bc)
{
  ggls<-list("with MRI stage + CRS + age + tumour.groupingII"=survfit(model11nostrat),"with CRS + age + tumour.groupingII"=survfit(model11a),"with MRI stage + age + tumour.groupingII"=survfit(model11bc))
  ggsc<-ggsurvplot_combine(ggls,xlab="OS Time (days)",title="Comparing the Primary Predictors Strength (MRI stage vs CRS): Cox (PH) Regression Model Comparisons by Survival Curves Estimate",legend.title = "Model by predictors", fin.df2)#, palette = "Dark2")
  ggsc$plot<-ggsc$plot+scale_colour_manual(values = c(rgb(1,0,0,1),rgb(0,0,1,0.5),rgb(0.9,0.65,0)))
  return(ggsc$plot)
}

plotIDa.predstrength.survCurvComp<-function(subpop){
  ggls<-list("with MRI stage + CRS + age + tumour.groupingII"=survfit(model11all,subset=subpop),"with CRS + age + tumour.groupingII"=survfit(model11alla,subset=subpop),"with MRI stage + age + tumour.groupingII"=survfit(model11allb,subset=subpop))
  ggsc<-ggsurvplot_combine(ggls,xlab="OS Time (days)",title="Comparing the Primary Predictors Strength (MRI stage vs CRS): Cox (PH) Regression Model Comparisons by Survival Curves Estimate",legend.title = "Model by predictors", fin.df2[subpop,])#, palette = "Dark2")
  ggsc$plot<-ggsc$plot+scale_colour_manual(values = c(rgb(1,0,0,1),rgb(0,0,1,0.5),rgb(0.9,0.65,0)))
  return(ggsc$plot)
}

subsetting<-function(){

summary(subset(model11nostrat,stage.MRI1=="IB"))
png("MRIstage&CRSinOS-subpopIB-subs.png",width=16,height=8,units='in',res=300)
grid.arrange(plotID.predstrength.survCurvComp(stage.MRI1=="IB"),bottom=paste("\nn = ",length(stageIBsubpop)," (subset: MRI stage = IB)"))
dev.off()

png("MRIstage&CRSinOS-subpopIB-subsid.png",width=16,height=8,units='in',res=300)
grid.arrange(plotID.predstrength.survCurvComp(stageIBsubpop),bottom=paste("\nn = ",length(stageIBsubpop)," (subset: MRI stage = IB)"))
dev.off()

png("MRIstage&CRSinOS-subpopIB-corr1.png",width=16,height=8,units='in',res=300)
gglss<-list("with MRI stage + CRS + age + tumour.groupingII"=survfit(model11nostrat,data=subset(fin.df2,stage.MRI1=="IB")),"with CRS + age + tumour.groupingII"=survfit(model11a,data=subset(fin.df2,stage.MRI1=="IB")),"with MRI stage + age + tumour.groupingII"=survfit(model11bc,data=subset(fin.df2,stage.MRI1=="IB")))
ggsc<-ggsurvplot_combine(gglss,xlab="OS Time (days)",title="Comparing the Primary Predictors Strength (MRI stage vs CRS): Cox (PH) Regression Model Comparisons by Survival Curves Estimate",legend.title = "Model by predictors", fin.df2)#, palette = "Dark2")
ggsc$plot<-ggsc$plot+scale_colour_manual(values = c(rgb(1,0,0,1),rgb(0,0,1,0.5),rgb(0.9,0.65,0)))
grid.arrange(ggsc$plot,bottom=paste("\nn = ",length(stageIBsubpop)," (subset: MRI stage = IB)"))
dev.off()

png("MRIstage&CRSinOS-subpopIB-corr2.png",width=16,height=8,units='in',res=300)
gglss<-list("with MRI stage + CRS + age + tumour.groupingII"=survfit(model11nostrat,data=subset(fin.df2,stage.MRI1=="IB")),"with CRS + age + tumour.groupingII"=survfit(model11a,data=subset(fin.df2,stage.MRI1=="IB")),"with MRI stage + age + tumour.groupingII"=survfit(model11bc,data=subset(fin.df2,stage.MRI1=="IB")))
ggsc<-ggsurvplot_combine(gglss,xlab="OS Time (days)",title="Comparing the Primary Predictors Strength (MRI stage vs CRS): Cox (PH) Regression Model Comparisons by Survival Curves Estimate",legend.title = "Model by predictors", data=subset(find.df2,stage.MRI1=="IB"))#, palette = "Dark2")
ggsc$plot<-ggsc$plot+scale_colour_manual(values = c(rgb(1,0,0,1),rgb(0,0,1,0.5),rgb(0.9,0.65,0)))
grid.arrange(ggsc$plot,bottom=paste("\nn = ",length(stageIBsubpop)," (subset: MRI stage = IB)"))
dev.off()

png("plot1.png",width=16,height=8,units='in',res=300)
ggsurvplot(survfit(coxph(ec.surv1 ~ stage.MRI1 + clinicrs + age1 + tumourv.groups2),data=fin.df2))
dev.off()
png("plot2.png",width=16,height=8,units='in',res=300)
ggsurvplot(survfit(coxph(ec.surv1 ~ stage.MRI1 + clinicrs + age1 + tumourv.groups2),data=subset(fin.df2,stage.MRI1=="IB")))
dev.off()

model11b<-coxph(ec.survb~crsb+ageb+tumourv.groupsb)
model11ab<-coxph(ec.survb~crsb+ageb+tumourv.groupsb)
model11bb<-coxph(ec.survb~ageb+tumourv.groupsb)
png("MRIstage&CRSinOS-subpopIB.png",width=16,height=8,units='in',res=300)
grid.arrange(plotM.predstrength.survCurvComp(model11b,model11ab,model11bb),bottom=paste("\nn = ",length(stageIBsubpop)," (subset: MRI stage = IB)"))
dev.off()

model11all<-coxph(ec.survall~stage.MRIall+crsall+ageall+tumourv.groupsall)
model11alla<-coxph(ec.survall~crsall+ageall+tumourv.groupsall)
model11allb<-coxph(ec.survall~stage.MRIall+ageall+tumourv.groupsall)
png("MRIstage&CRSinOS-allpop.png",width=16,height=8,units='in',res=300)
grid.arrange(plotM.predstrength.survCurvComp(model11all,model11alla,model11allb),bottom=paste("\nn = ",length(allstudypop)," (our study population)"))
dev.off()

#grid.arrange(predstrength.survCurvComp(1:3),bottom=paste("\nn = ",50," (subset: MRI stage = IB)"))
model11n4<-coxph(ec.survn4~stage.MRIn4+crsn4+agen4+tumourv.groupsn4)
model11an4<-coxph(ec.survn4~crsn4+agen4+tumourv.groupsn4)
model11nb4<-coxph(ec.survn4~stage.MRIn4+agen4+tumourv.groupsn4)
png("MRIstage&CRSinOS-subpopWoutIV.png",width=16,height=8,units='in',res=300)
grid.arrange(plotM.predstrength.survCurvComp(model11n4,model11an4,model11nb4),bottom=paste("\nn = ",length(stagenIVsubpop)," (subset: MRI stage = all except for IV)"))
dev.off()

png("MRIstage&CRSinOS-allpop-subpopWoutIV.png",width=16,height=8,units='in',res=300)
grid.arrange(plotIDa.predstrength.survCurvComp(stagenIVsubpop),bottom=paste("\nn = ",length(stagenIVsubpop)," (subset: MRI stage = all except for IV)"))
dev.off()

png("MRIstage&CRSinOS-allpop-subpopIB.png",width=16,height=8,units='in',res=300)
grid.arrange(plotIDa.predstrength.survCurvComp(which(stage.MRIall=="IB")),bottom=paste("\nn = ",length(stageIBsubpop)," (subset: MRI stage = IB)"))
dev.off()
}


# we decide on model 11
ultimate.comparison<-function(){
#in stage I patients
#models 6a and 6b are not good, even though model 6 is
#model6a<-coxph(ec.surv1~clinicrs+age.groups1+strata(tumourv.groups2))
stage.MRIs<-strata(stage.MRI1)
model11a<-coxph(ec.surv1 ~ clinicrs + age1 + 
                  tumourv.groups2)
model11b<-coxph(ec.surv1 ~ stage.MRIs + age1 + 
                  tumourv.groups2)
# model10a<-coxph(ec.surv1 ~ clinicrs + ages + tumourv.groups2)
# model10b<-coxph(ec.surv1 ~ stage.MRI1 + ages + tumourv.groups2)
# # model6b<-coxph(ec.surv1~stage.MRI1+age.groups1+strata(tumourv.groups2))
fit4.df<-data.frame(ec.surv1[,1],ec.surv1[,2],clinicrs,stage.MRI1,age1,tumour.vol1,tumourv.groups2)
fit4 <- list(complete.with.both = survfit(model11), without.MRI.stage = survfit(model11a),without.CRS = survfit(model11b), MRI.stage.better.fit=survfit(stratifmodel))
png("MRIstage&CRSinOS-model11-stratif.png",width=16,height=8,units='in',res=300)
ggsurvplot_combine(fit4,xlab="OS Time (days)",title="Comparison between estimated survival curves for different predictors of OS: MRI stage vs CRS",subtitle="(with the complete model stratified by MRI Stage)",legend.title = "Model by primary predictors", fit4.df, palette = "Dark2")
dev.off()
model11all<-coxph(ec.surv1 ~ stage.MRIs + clinicrs + age1 + 
        tumourv.groups2)
fit4all <- list(complete.with.both = survfit(model11all), without.MRI.stage = survfit(model11a),without.CRS = survfit(model11b), MRI.stage.better.fit=survfit(stratifmodel))
png("MRIstage&CRSinOS-model11.png",width=16,height=8,units='in',res=300)
ggsurvplot_combine(fit4all,xlab="OS Time (days)",title="Comparison between estimated survival curves for different predictors of OS: MRI stage vs CRS",legend.title = "Model by primary predictors",legend.labs =c("with MRI stage and CRS","without MRI stage","without CRS","with MRI stage and histo grade"), fit4.df)#, palette = "Dark2")
dev.off()

stratmodelIBa<-coxph(ec.survb ~ crsb + ageb + tumourv.groupsb)
stratmodelIBb<-coxph(ec.survb ~ strata(stage.MRIb) + ageb + tumourv.groupsb)
fit1b<-list(complete.with.both = survfit(stratmodelIB), without.MRI.stage = survfit(stratmodelIBa),without.CRS = survfit(stratmodelIBb))
fitb.df<-data.frame(ec.survb[,1],ec.survb[,2],crsb,stage.MRIb,ageb,tumour.volb,tumourv.groupsb)
png("IBMRIstage&CRSinOS-model11.png",width=16,height=8,units='in',res=300)
ggsurvplot_combine(fit1b,xlab="OS Time (days)",title="Comparison between estimated survival curves for different predictors of OS: MRI stage vs CRS",subtitle="Model 11 ajusted for stage IB subpopulation",legend.title = "Model by primary predictors",legend.labs =c("with MRI stage and CRS","without MRI stage","without CRS"), fitb.df)
dev.off()

# sall<-strata(stage.MRIall)
# stratamodelall<-coxph(ec.survall ~ strata(stage.MRIall) + crsall + ageall + tumourv.groupsall)
# stratmodelalla<-coxph(ec.survall ~ crsall + ageall + tumourv.groupsall)
# stratmodelallb<-coxph(ec.survall ~ strata(stage.MRIall) + ageall + tumourv.groupsall)
# #fitall<-list(complete.with.both = survfit(stratmodelall), 
# fitall<-list(without.MRI.stage = survfit(stratmodelalla),without.CRS = survfit(stratmodelallb))#,study.better.fit=survfit(model13))
# fitall.df<-data.frame(ec.survall[,1],ec.survall[,2],crsall,stage.MRIall,ageall,tumour.volall,tumourv.groupsall)
# png("MRIstage&CRSinOS-model11-studypop.png",width=16,height=8,units='in',res=300)
# ggsurvplot_combine(fitall,xlab="OS Time (days)",title="Comparison between estimated survival curves for different predictors of OS: MRI stage vs CRS",subtitle="Model 11 ajusted for our study population",legend.title = "Model by primary predictors",data=fitall.df,legend.labs =c("with MRI stage and CRS","without MRI stage","without CRS"))#, fitall.df)
# dev.off()

plot_summs(model11,model13)
}


# full models (6 or 11) vs Kaplan-Meier with only one predictor

modelbyMRI<-function(){
#plot
#surv.median.line = "hv"
fit<-survfit(ec.surv1~stage.MRI1+age1+tumourvolstrat+grade.histo1)
fit1<-survfit(ec.surv1~stage.MRI1+age1+tumourvolstrat+grade.histo1,subset = stage.MRI1)
surv.df<-data.frame(ec.surv1,stage.MRI1,age1,tumourvolstrat,grade.histo1)
# ggsurvs<-ggsurvplot(fit,data=surv.df, fun = "event", conf.int = TRUE,
#            ggtheme = theme_bw())#,facet.by ="stage.MRI1")
# 
# ggsurvs$plot +theme_bw() +theme (legend.position = "right")+facet_grid(age.groups1~stage.MRI1)

ggsurvplot(survfit(ec.surv1~stage.MRI1),data=surv.df, fun = "event", conf.int = TRUE,
           ggtheme = theme_bw())
# df2<-data.frame(stage.MRI1,clinicrs)
# ggsurvplot(fit,df2, conf.int = TRUE,
#            ggtheme = theme_bw(),facet.by = "clinicrs")
clinicrs<-corrected.data$Risk.score[early.st.pop]

png("SurvivalCurvesByMRIstage-FULLvs1pred.png",width=16,height=8,units='in',res=300)
ggs<-ggsurvplot(survfit(ec.surv1~stage.MRI1),data=surv.df,xlab="Time (days)",title="Kaplan-Meier estimate with MRI Stage unadjusted",risk.table=TRUE,legend.title="MRI Stage",legend.labs=levels(stage.MRI1),
           ggtheme = theme_bw())
gga<-ggadjustedcurves(stratmodel3,variable = "stage.MRI1",data = surv.df1,reference = NULL,method = "average",fun = NULL,title="Cox (PH) model estimate with MRI stage adjusted for\nage + strata(tumour.groupingII) + tumour grade",xlab="Time (days)",font.title=c(14,"plain","black"),font.x=c(11,"plain","black"),font.y=c(11,"plain","black"),font.legend=c(9,"plain","black"))
gga<-gga+ theme(axis.text.x=element_text(size=10,color="grey40"),axis.text.y=element_text(size=10,color="grey40"))
#gga<-gga+labs(title="Based on Cox Model and Separated by MRI Stage",x="Time (days)")
grid.arrange(ggs$plot,gga,ggs$table,
             layout_matrix = rbind(c(1,2),c(1,2),c(1,2),c(3,NA)),top="Survival Curves\n")#c(1, 2, NA),c(3, 3, 4)
dev.off()

library(ggpubr)
library(ggprism)
strattumourvg2<-strata(tumourv.groups2)
fit6<-survfit(ec.surv1 ~ stage.MRI1 + clinicrs + age.groups1 + strata(tumourv.groups2),data=surv.df2)
stat.test<-surv_pvalue(fit6,method="survdiff")
#p < 0.0001
survdiff(ec.surv1 ~ stage.MRI1 + clinicrs + age.groups1 + strattumourvg2,data=surv.df2)
#p= 2e-11
fit11<-survfit(ec.surv1 ~ strata(stage.MRI1) + clinicrs + age1 + 
                 tumourv.groups2,data=fit4.df)
stat.test1<-surv_pvalue(fit11,method="survdiff")
#p < 0.0001
stat.test1$pval
survdiff(ec.surv1 ~ strata(stage.MRI1) + clinicrs + age1 + 
           tumourv.groups2,data=fit4.df)
#p= <2e-16 
#surv.dfg<-data.frame(group1=stage.MRI1[stage.MRI1=="IA"],group2=stage.MRI1[stage.MRI1=="IB"],label=stat.test,y.position=c(0.25,500))
png("SurvivalCurvesByMRIstage-FULLvs1pred-commonModelPval-model11-sprob.png",width=16,height=8,units='in',res=300)
title<-grid.text("The Overall Survival Grouped by the MRI Stage",gp=gpar(fontsize=20,fontface="bold",fontfamily="Calibri"))#,fontfamily="Calibri"))
subtitle<-grid.text("Survival Curves\n",gp=gpar(fontsize=16,fontfamily="Arial"))
margin <- unit(0.5, "line")
ggs<-ggsurvplot(survfit(ec.surv1~stage.MRI1),data=surv.df,xlab="OS Time (days)",title="Kaplan-Meier estimate with MRI Stage unadjusted",risk.table=TRUE,legend.title="MRI Stage",legend.labs=levels(stage.MRI1))
                #pval = TRUE)#ggtheme = theme_bw(),
p<-surv_pvalue(
  fit=survfit(ec.surv1~stage.MRI1),
  data = fit4.df,
  method = "survdiff",
  test.for.trend = FALSE,
  combine = TRUE
)
pval<-p$pval
ggs$plot<-ggs$plot+annotate("text", x = 750, y=0.2, label = paste("p = ",pval," for log-rank test"),size=5)
model11nostrat<-coxph(formula = ec.surv1 ~ stage.MRI1 + clinicrs + age1 + 
                        tumourv.groups2)
p<-surv_pvalue(
  fit=survfit(model11nostrat$formula,data=fit4.df),
  data = fit4.df,
  method = "survdiff",
  test.for.trend = FALSE,
  combine = TRUE
)
pval<-p$pval
ggatab<-ggsurvtable(survfit(model11nostrat),data=fit4.df, timeby = 1000,xlab="OS Time (days)")
gga<-ggadjustedcurves(model11nostrat,variable = "stage.MRI1",data = fit4.df,method = "average",title="Survival curves predicted by the Cox (PH) model with MRI stage adjusted\nfor CRS + age + tumour.groupingII",xlab="OS Time (days)",ylab="Survival Probability",legend="top",legend.title="MRI stage")#,font.title=c(14,"plain","black"),font.x=c(11,"plain","black"),font.y=c(11,"plain","black"),font.legend=c(9,"plain","black"))
gga<-gga+ annotate("text", x = 750, y=0.4, label = paste("p = ",pval," for log-rank test"),size=5)#+stat_pvalue_manual(surv.df2,y.position = "3",xmin=stage.MRI1[stage.MRI1=="IA"],xmax=stage.MRI1[stage.MRI1=="IB"],label="p")#add_pvalue(surv.dfg,label="p = {label}")#,y.position = c(0.25,500))
#theme(axis.text.x=element_text(size=10,color="grey40"),axis.text.y=element_text(size=10,color="grey40"))
#gga<-gga+labs(title="Based on Cox Model and Separated by MRI Stage",x="Time (days)")
gga<-gga+annotate("text", x = 750, y=0.1, label = paste("HR = 0.863 for MRI stage (IB relative to IA)\n95% CI: 0.4840-1.539"),size=5)
#gga<-gga+ggatab$risk.table
btm<-grid.text("\nn = 319 (subset: MRI stage = I; IA-IB = unbalanced groups)",gp=gpar(fontsize=16,fontfamily="Arial"))
gra<-grid.arrange(ggs$plot,gga,ggs$table,ggatab$risk.table,layout_matrix = rbind(c(1,2),c(1,2),c(1,2),c(3,4)))
grid.arrange(title,subtitle,gra,heights = unit.c(grobHeight(title) + 1.2*margin, 
                                              grobHeight(subtitle) + margin, 
                                              unit(1,"null")), bottom=btm)
dev.off()

png("CoxModelGrByMRIst&histo.png",width=16,height=8,units='in',res=300)
ggm<-ggadjustedcurves(model11nostrat,variable = "stage.MRI1",data = MRIhisto.df,method = "average",title="grouped by MRI stage",xlab="OS Time (days)",ylab="Survival Probability",legend="top",legend.title="MRI stage")#,font.title=c(14,"plain","black"),font.x=c(11,"plain","black"),font.y=c(11,"plain","black"),font.legend=c(9,"plain","black"))
ggh<-ggadjustedcurves(model11nostrat,variable = "stage.histo1",data = MRIhisto.df,method = "average",title="grouped by histological stage",xlab="OS Time (days)",ylab="Survival Probability",legend="top",legend.title="histo stage")

title<-grid.text("Survival Curves predicted by the Cox model with MRI stage + CRS + age + tumour.groupingII",gp=gpar(fontsize=20,fontface="bold",fontfamily="Calibri"))
ggt<-grid.arrange(ggm,ggh,ncol=2)
grid.arrange(title,ggt,heights = unit.c(grobHeight(title) + 1.2*margin,unit(1,"null")),bottom="\nn = 319 (subset: MRI stage = I; unbalanced groups for stages)")
dev.off()

#arrange_ggsurvplots(list("by MRI stage"=ggs1,"by histological stage"=ggs2),data = MRIhisto.df,title="Survival curves predicted by the Cox (PH) model with MRI stage adjusted\nfor CRS + age + tumour.groupingII",xlab="OS Time (days)",ylab="Survival Probability",legend="top",legend.title="Survival Curves Grouped")
png("SurvivalCurvesByMRI&histoStage-FULLvs1pred-commonModelPval-model11-sprob.png",width=14,height=10,units='in',res=300)
gg<-ggsurvplot_combine(list("MRI stage"=survfit(ec.surv1~stage.MRI1),"histological stage"=survfit(ec.surv1~stage.histo1)),data = MRIhisto.df,risk.table = TRUE,title="Two Kaplan-Meier Models plotted together\n(one with MRI stage alone, the other with histological stage alone)",xlab="OS Time (days)",ylab="Survival Probability",legend="top",legend.title="Survival Curves Grouped by",legend.labs=c("MRI stage:\nIA","IB ","histo stage:\nIA","IB","II","IIIA" ,"IIIC1","IIIC2","IV"))

#gg + labs(stage.MRI1="MRI stage", stage.histo1 ="histological stage")

#theme(axis.text.x=element_text(size=10,color="grey40"),axis.text.y=element_text(size=10,color="grey40"))
#gga<-gga+labs(title="Based on Cox Model and Separated by MRI Stage",x="Time (days)")
#gga<-gga+annotate("text", x = 750, y=0.1, label = paste("HR = 0.863 for MRI stage (IB relative to IA)\n95% CI: 0.4840-1.539"),size=5)
#gga<-gga+ggatab$risk.table
btm<-grid.text("\nn = 319 (subset: MRI stage = I; unbalanced groups for stages)",gp=gpar(fontsize=16,fontfamily="Arial"))
#gra<-grid.arrange(ggs$plot,gga,ggs$table,ggatab$risk.table,layout_matrix = rbind(c(1,2),c(1,2),c(1,2),c(3,4)))
grid.arrange(gg$plot,gg$table,btm,layout_matrix = rbind(c(1,1),c(1,1),c(1,1),c(1,1),c(1,1),c(1,1),c(2,2),c(2,2),c(2,2),c(3,3)))
#grid.arrange(title,subtitle,ggr,heights = unit.c(grobHeight(title) + 1.2*margin, 
                                                 # grobHeight(subtitle) + margin, 
                                                 # unit(1,"null")), bottom=btm)
dev.off()

png("SurvivalCurvesByMRIstage-FULLvs1pred-commonModelPval-model6.png",width=16,height=8,units='in',res=300)
ggs<-ggsurvplot(survfit(ec.surv1~stage.MRI1),data=cbind(surv.df2,ec.surv1[,2]),xlab="OS Time (days)",title="Kaplan-Meier estimate with MRI Stage unadjusted",risk.table=TRUE,legend.title="MRI Stage",legend.labs=levels(stage.MRI1),
                pval = TRUE)
p<-surv_pvalue(
  fit=surv_fit(model6$formula,data=surv.df),
  data = surv.df,
  method = "survdiff",
  test.for.trend = FALSE,
  combine = TRUE
)
pval<-p$pval
gga<-ggadjustedcurves(model6,variable = "stage.MRI1",data=cbind(surv.df2,ec.surv1[,2]),method = "average",title="Cox (PH) model estimate with MRI stage adjusted for\nCRS + age + tumour.groupingII",xlab="OS Time (days)",ylab="Survival Probability",legend="top",legend.title="MRI stage")#title="Based on Cox Model and Separated by MRI Stage",xlab="Time (days)",font.title=c(14,"plain","black"),font.x=c(11,"plain","black"),font.y=c(11,"plain","black"),font.legend=c(9,"plain","black"))
gga<-gga+ annotate("text", x = 500, y=0.2, label = paste("p = ",pval),size=5)#+stat_pvalue_manual(surv.df2,y.position = "3",xmin=stage.MRI1[stage.MRI1=="IA"],xmax=stage.MRI1[stage.MRI1=="IB"],label="p")#add_pvalue(surv.dfg,label="p = {label}")#,y.position = c(0.25,500))
#gga<-gga+labs(title="Based on Cox Model and Separated by MRI Stage",x="Time (days)")
gga<-gga+annotate("text", x = 500, y=0.1, label = paste("HR = 1.011 for MRI stage"),size=5)
btm<-grid.text("\nn = 319 (subset: MRI stage = I; IA-IB = unbalanced groups)",gp=gpar(fontsize=16,fontfamily="Arial"))
gra<-grid.arrange(ggs$plot,gga,ggs$table,ggatab$risk.table,layout_matrix = rbind(c(1,2),c(1,2),c(1,2),c(3,4)))
grid.arrange(ggs$plot,gga,ggs$table,
             layout_matrix = rbind(c(1,2),c(1,2),c(1,2),c(3,NA)),top="Survival Curves\n")#c(1, 2, NA),c(3, 3, 4)
dev.off()

png("SurvivalCurvesByMRIstage-FULLvs1pred-commonModel-model11.png",width=16,height=8,units='in',res=300)
ggs<-ggsurvplot(survfit(ec.surv1~stage.MRI1),data=surv.df,xlab="Time (days)",title="KM Model with Only One Predictor (MRI Stage)",risk.table=TRUE,legend.title="MRI Stage",legend.labs=levels(stage.MRI1),
                ggtheme = theme_bw())
gga<-ggadjustedcurves(model11,variable = "stage.MRI1",data = fit4.df,reference = NULL,method = "average",fun = NULL,title="Based on Cox Model and Separated by MRI Stage",xlab="Time (days)",font.title=c(14,"plain","black"),font.x=c(11,"plain","black"),font.y=c(11,"plain","black"),font.legend=c(9,"plain","black"))
gga<-gga+ theme(axis.text.x=element_text(size=10,color="grey40"),axis.text.y=element_text(size=10,color="grey40"))
grid.arrange(ggs$plot,gga,ggs$table,
             layout_matrix = rbind(c(1,2),c(1,2),c(1,2),c(3,NA)),top="Survival Curves\n")#c(1, 2, NA),c(3, 3, 4)
dev.off()
}

modelbyCRS<-function()
{
png("SurvivalCurvesByCRS-FULLvs1pred-commonModel-model11.png",width=16,height=8,units='in',res=300)
clinicrs1<-clinicrs
clinicrs<-factor(clinicrs1,levels=c("low","intermediate","high","advanced"))
levels(clinicrs)<-rev(c("advanced","high","intermediate","low"))
ggs<-ggsurvplot(
  survfit(ec.surv1 ~ clinicrs),
  title = "The Overall Survival Based on Clinical Risk Score",
  #subtitle = "Stratification by Risk Score",
  data = data,
  legend.title = "CR Score",
  risk.table = TRUE,
  #pval = TRUE,
  conf.int = FALSE, # show confidence intervals for 
  # point estimates of survival curves.
  #xlim = c(0,5000),
  xlab = "OS Time (days)",
  break.time.by = 500,
  #ggtheme = theme_light(),
  #brewer.pal(4,"Set2")[3]
  palette = rev(c("red","orange","deepskyblue","forestgreen")),
  risk.table.title="Strata (by CR Score) - size every 500 days",
  #risk.table.ticks.col = TRUE,
  #risk.table.y.text.col = T,
  #risk.table.y.text = FALSE,
  ggtheme = theme_bw(),
  legend.labs = levels(clinicrs)#c("advanced","high","intermediate","low")
  #legend.labs = c("risk.score=advanced", "risk.score=high", "risk.score=intermediate", "risk.score=low")
)
#model6<-coxph(ec.surv1 ~ stage.MRI1 + clinicrs + age.groups1 + 
                #strata(tumourv.groups2))
gga<-ggadjustedcurves(model11,variable = "clinicrs",data = fit4.df,reference = NULL,method = "average",fun = NULL,title="Based on Cox Model and Separated by Clinical RS",xlab="Time (days)",font.title=c(14,"plain","black"),font.x=c(11,"plain","black"),font.y=c(11,"plain","black"),font.legend=c(9,"plain","black"),palette = rev(c("red","orange","deepskyblue","forestgreen")))
gga<-gga+ theme(axis.text.x=element_text(size=10,color="grey40"),axis.text.y=element_text(size=10,color="grey40"))#+ annotate("text", x = 500, y=0.2, label = paste("p = ",stat.test$pval),size=5,fontface='bold')
grid.arrange(ggs$plot,gga,ggs$table,
             layout_matrix = rbind(c(1,2),c(1,2),c(1,2),c(3,NA)),top="Survival Curves\n")#c(1, 2, NA),c(3, 3, 4)
dev.off()

png("SurvivalCurvesByCRS-FULLvs1pred-commonModel-model6.png",width=16,height=8,units='in',res=300)
ggs<-ggsurvplot(
  survfit(ec.surv1 ~ clinicrs),
  title = "The Overall Survival Based on Clinical Risk Score",
  #subtitle = "Stratification by Risk Score",
  data = data,
  legend.title = "CR Score",
  risk.table = TRUE,
  #pval = TRUE,
  conf.int = FALSE, # show confidence intervals for 
  # point estimates of survival curves.
  #xlim = c(0,5000),
  xlab = "OS Time (days)",
  break.time.by = 500,
  #ggtheme = theme_light(),
  #brewer.pal(4,"Set2")[3]
  palette = rev(c("red","orange","deepskyblue","forestgreen")),
  risk.table.title="Strata (by CR Score) - size every 500 days",
  #risk.table.ticks.col = TRUE,
  #risk.table.y.text.col = T,
  #risk.table.y.text = FALSE,
  ggtheme = theme_bw(),
  legend.labs = levels(clinicrs)#c("advanced","high","intermediate","low")
  #legend.labs = c("risk.score=advanced", "risk.score=high", "risk.score=intermediate", "risk.score=low")
)
#model6<-coxph(ec.surv1 ~ stage.MRI1 + clinicrs + age.groups1 + 
#strata(tumourv.groups2))
gga<-ggadjustedcurves(model6,variable = "clinicrs",data = fit4.df,reference = NULL,method = "average",fun = NULL,title="Based on Cox Model and Separated by Clinical RS",xlab="Time (days)",font.title=c(14,"plain","black"),font.x=c(11,"plain","black"),font.y=c(11,"plain","black"),font.legend=c(9,"plain","black"),palette = rev(c("red","orange","deepskyblue","forestgreen")))
gga<-gga+ theme(axis.text.x=element_text(size=10,color="grey40"),axis.text.y=element_text(size=10,color="grey40"))#+ annotate("text", x = 500, y=0.2, label = paste("p = ",stat.test$pval),size=5,fontface='bold')
grid.arrange(ggs$plot,gga,ggs$table,
             layout_matrix = rbind(c(1,2),c(1,2),c(1,2),c(3,NA)),top="Survival Curves\n")

dev.off()

png("SurvivalCurvesByCRS-FULLvs1pred-commonModelPval-model11.png",width=16,height=8,units='in',res=300)
ggs<-ggsurvplot(
  survfit(ec.surv1 ~ clinicrs),
  title = "Kaplan-Meier estimate with CRS unadjusted",
  #subtitle = "Stratification by Risk Score",
  data = data,
  legend.title = "CR Score",
  risk.table = TRUE,
  #pval = TRUE,
  conf.int = FALSE, # show confidence intervals for 
  # point estimates of survival curves.
  #xlim = c(0,5000),
  xlab = "OS Time (days)",
  break.time.by = 500,
  #ggtheme = theme_light(),
  #brewer.pal(4,"Set2")[3]
  palette = rev(c("red","orange","deepskyblue","forestgreen")),
  risk.table.title="Strata (by CR Score) - size every 500 days",
  #risk.table.ticks.col = TRUE,
  #risk.table.y.text.col = T,
  #risk.table.y.text = FALSE,
  #ggtheme = theme_bw(),
  legend.labs = levels(clinicrs)#c("advanced","high","intermediate","low")
  #legend.labs = c("risk.score=advanced", "risk.score=high", "risk.score=intermediate", "risk.score=low")
)
#model6<-coxph(ec.surv1 ~ stage.MRI1 + clinicrs + age.groups1 + 
#strata(tumourv.groups2))
ggatab<-ggsurvtable(survfit(model11nostrat),data=fit4.df, timeby = 1000,xlab="OS Time (days)")
gga<-ggadjustedcurves(model11nostrat,variable = "clinicrs",data = fit4.df,method = "average",palette = rev(c("red","orange","deepskyblue","forestgreen")),title="Survival curves predicted by the Cox (PH) model with CRS adjusted\nfor MRI stage + age + tumour.groupingII",xlab="OS Time (days)",ylab="Survival Probability",legend="top",legend.title="CRS")
#font.title=c(14,"plain","black"),font.x=c(11,"plain","black"),font.y=c(11,"plain","black"),font.legend=c(9,"plain","black")
p<-surv_pvalue(
  fit=surv_fit(model11nostrat$formula,data=fit4.df),
  data = fit4.df,
  method = "survdiff",
  test.for.trend = FALSE,
  combine = TRUE
)
pval<-p$pval
gga<-gga+ annotate("text", x = 750, y=0.2, label = paste("p = ",pval," for log-rank test"),size=5)#fontface='bold')
#gga<-gga+annotate("text", x = 500, y=0.1, label = paste("HR = 0.863 for MRI stage"),size=5)
#theme(axis.text.x=element_text(size=10,color="grey40"),axis.text.y=element_text(size=10,color="grey40"))+
btm<-grid.text("\nn = 319 (subset: MRI stage = I; IA-IB = unbalanced groups)",gp=gpar(fontsize=16,fontfamily="Arial"))
title<-grid.text("The Overall Survival Grouped by the Clinical Risk Score",gp=gpar(fontsize=20,fontface="bold",fontfamily="Calibri"))#,fontfamily="Calibri"))
subtitle<-grid.text("Survival Curves\n",gp=gpar(fontsize=16,fontfamily="Arial"))
margin <- unit(0.5, "line")
gra<-grid.arrange(ggs$plot,gga,ggs$table,ggatab$risk.table,layout_matrix = rbind(c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(3,4),c(3,4)))
grid.arrange(title,subtitle,gra,heights = unit.c(grobHeight(title) + 1.2*margin, 
                                                 grobHeight(subtitle) + margin, 
                                                 unit(1,"null")), bottom=btm)#c(1, 2, NA),c(3, 3, 4)
dev.off()

png("SurvivalCurvesByCRS-FULLvs1pred-commonModelPval-model6.png",width=16,height=8,units='in',res=300)
ggs<-ggsurvplot(
  survfit(ec.surv1 ~ clinicrs),
  title = "The Overall Survival Based on Clinical Risk Score",
  #subtitle = "Stratification by Risk Score",
  data = data,
  legend.title = "CR Score",
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE, # show confidence intervals for 
  # point estimates of survival curves.
  #xlim = c(0,5000),
  xlab = "OS Time (days)",
  break.time.by = 500,
  #ggtheme = theme_light(),
  #brewer.pal(4,"Set2")[3]
  palette = rev(c("red","orange","deepskyblue","forestgreen")),
  risk.table.title="Strata (by CR Score) - size every 500 days",
  #risk.table.ticks.col = TRUE,
  #risk.table.y.text.col = T,
  #risk.table.y.text = FALSE,
  ggtheme = theme_bw(),
  legend.labs = levels(clinicrs)#c("advanced","high","intermediate","low")
  #legend.labs = c("risk.score=advanced", "risk.score=high", "risk.score=intermediate", "risk.score=low")
)
#model6<-coxph(ec.surv1 ~ stage.MRI1 + clinicrs + age.groups1 + 
#strata(tumourv.groups2))
gga<-ggadjustedcurves(model6,variable = "clinicrs",data = fit4.df,reference = NULL,method = "average",fun = NULL,title="Based on Cox Model and Separated by Clinical RS",xlab="Time (days)",font.title=c(14,"plain","black"),font.x=c(11,"plain","black"),font.y=c(11,"plain","black"),font.legend=c(9,"plain","black"),palette = rev(c("red","orange","deepskyblue","forestgreen")))
gga<-gga+ theme(axis.text.x=element_text(size=10,color="grey40"),axis.text.y=element_text(size=10,color="grey40"))+ annotate("text", x = 500, y=0.2, label = paste("p = ",stat.test$pval),size=5,fontface='bold')
grid.arrange(ggs$plot,gga,ggs$table,
             layout_matrix = rbind(c(1,2),c(1,2),c(1,2),c(3,NA)),top="Survival Curves\n")

dev.off()

#to plot
surv.df1<-data.frame(ec.surv1,stage.MRI1,age1,tumourv.groups2,grade.histo1)
ggadjustedcurves(stratmodel3,variable = "stage.MRI1",data = surv.df1,reference = NULL,method = "conditional",fun = NULL)
ggadjustedcurves(stratmodel3,variable = "stage.MRI1",data = surv.df1,reference = NULL,method = "average",fun = NULL)
#ggcoxadjustedcurves(stratmodel3,data = surv.df1,variable = surv.df1[,1])
plot(survfit(ec.surv1 ~ 1, data = surv.df1))
#ggsurvplot(fit1, data = surv.df1)

llrs<-list()
llrs["low"]<-length(which(clinicrs=="low"))
llrs["intermediate"]<-length(which(clinicrs=="intermediate"))
llrs["high"]<-length(which(clinicrs=="high"))
llrs["advanced"]<-length(which(clinicrs=="advanced"))
lrs<-paste0(levels(clinicrs),"\n","n = ",llrs[levels(clinicrs)])
png("CRS&OS-correlation.png",width=16,height=8,units='in',res=300)
var(ec.surv1[,1])
hist(ec.surv1[,1])
ecsurv1.time<-ec.surv1[,1]
#wilcox.test(ecsurv1.time~clinicrs,paired=TRUE);cannot accommodate more than two groups
ggplot(surv.df1,aes(x = clinicrs, y = ecsurv1.time))+   geom_boxplot(size = .75, outlier.shape = NA) +stat_compare_means(paired=TRUE,label.x=3.95,label.y=3000)+geom_jitter(aes(color=clinicrs),position=position_jitter(0.2))+  labs(x="Clinical Risk Score",y="Time (days)",title="Relationship between CRS and Survival Time (n = 319, subset: MRI stage = I)\n")+theme(legend.position = "none",plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+scale_x_discrete(labels=lrs)
#ggp2<-ggplot(surv.df1,aes(x = clinicrs, y = ec.surv1[,1], color = clinicrs))+   geom_jitter(position=position_jitter(0.2))+ labs(x="Clinical Risk Score",y="Time (days)")+theme(legend.position = "none",plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
#grid.arrange(ggp1,ggp2,ncol=2, top = "Relationship between CRS and Survival Time (n = 319, subset: MRI stage = I)\n")
dev.off()
png("MRIstage&OS-correlation.png",width=16,height=8,units='in',res=300)
#comparisons=list(stage.MRI1,ecsurv1.time)
library(varhandle)
stage.MRI1unf<-unfactor(stage.MRI1)
llst<-list()
llst["IA"]<-length(which(stage.MRI1=="IA"))
llst["IB"]<-length(which(stage.MRI1=="IB"))
# llst["II"]<-length(which(stage.MRI1=="II"))
# llst["IIIA"]<-length(which(stage.MRI1=="IIIA"))
# llst["IIIB"]<-length(which(stage.MRI1=="IIIB"))
# llst["IIIC1"]<-length(which(stage.MRI1=="IIIC1"))
# llst["IIIC2"]<-length(which(stage.MRI1=="IIIC2"))
# llst["IV"]<-length(which(stage.MRI=="IV"))
lst<-paste0(levels(stage.MRI),"\n","n = ",llst[levels(stage.MRI1)])
ggplot(surv.df1,aes(x = stage.MRI1unf, y = ecsurv1.time))+ geom_boxplot(size = .75, outlier.shape = NA) +stat_compare_means(method="kruskal.test",paired=TRUE,label.x=2.27,label.y=3000)+geom_jitter(aes(color=stage.MRI1unf),position=position_jitter(0.2)) + labs(x="MRI Stage",y="Time (days)",title = "Relationship between MRI Stage and Survival Time (n = 319, subset: MRI stage = I)\n")+theme(legend.position = "none",plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+scale_x_discrete(labels=lst)
#ggp2<-ggplot(surv.df1,aes(x = stage.MRI1, y = ec.surv1[,1], color = stage.MRI1))+   geom_jitter(position=position_jitter(0.2))+ labs(x="MRI Stage",y="Time (days)")+theme(legend.position = "none",plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
#grid.arrange(ggp1,ggp2, ncol = 2, top = "Relationship between MRI Stage and Survival Time (n = 319, subset: MRI stage = I)\n")
dev.off()
}
