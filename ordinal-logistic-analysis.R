logistic.regression<-function()
{
  table(corrected.data$TumourVolume.mm.3.[allstudypop])
  stageIpop<-intersect(which(corrected.data$X.2 %in% c("IA","IB")),allstudypop)
  volpop<-setdiff(allstudypop,which(is.na(corrected.data$TumourVolume.mm.3.)))
  #tumour volume is the exposure and it is continuous,
  #while MRI stage is the outcome and it is nominal
  stage.MRI<-factor(corrected.data$X.2[volpop])
  stage.MRI.scale<-as.numeric(stage.MRI)
  hist(stage.MRI.scale)
  sd(stage.MRI.scale)
  min(stage.MRI.scale)
  max(stage.MRI.scale)
  stage.MRI.rescale<-(stage.MRI.scale-min(stage.MRI.scale))/(max(stage.MRI.scale)-min(stage.MRI.scale))
  sd(stage.MRI.rescale)
  summary(stage.MRI.rescale)
  stage.MRI.rescalef<-cut(stage.MRI.rescale,8)
  table(stage.MRI.rescalef)
  table(stage.MRI)
  
  tumour.vol<-corrected.data$TumourVolume.mm.3.[volpop]
  hist(tumour.vol)
  tumour.vol.rescale<-(tumour.vol-min(tumour.vol))/(max(tumour.vol)-min(tumour.vol))
  sd(tumour.vol.rescale)
  summary(tumour.vol.rescale)
  length(stage.MRI.rescale)
  length(tumour.vol.rescale)
  
  tumour.volm3<-tumour.vol/(10^9)
  tumour.volm3.rescale<-(tumour.volm3-min(tumour.volm3))/(max(tumour.volm3)-min(tumour.volm3))
  sd(tumour.volm3.rescale)#same value when scaled (mm3 vs m3)
  
  #ggpairs(array.of.predictors, title = "Correlation Plot between each (predictor) variable")
  library(MASS)
  #polr is a special glm
  logistic.model<-polr(stage.MRI.rescalef~tumour.vol.rescale,Hess=TRUE,data=corrected.data)
  #since I intend to call summary on the fit, I should tell 
  #the function to return the Hessian (the observed information matrix) 
  summary(logistic.model)
  #the t-value is the regression coefficient to its standard error
  #intercepts indicate where the latent variable is cut to make the three groups that we observe in our data
  #intercepts are not used in the interpretation of the results,
  #whereas residual deviance and AIC are useful for interpretation of the results (model comparison)
  #residual deviance is equal to -2 * Log Likelihood of the model
  ctable <- coef(summary(logistic.model))
  #table of the regression coefficients
  p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
  ctable <- cbind(ctable, "p value" = p)
  #advanced stages have the largest absolute values of the coefficients;
  #all have p<0.05, meaning that they are statistically significant
  #we can check this with confindence intervals as well
  ci <- confint(logistic.model)
  confint.default(logistic.model) # CIs assuming normality
  #the conf.int does not contain 0, so the parameter estimate is stat.sign.
  #for a unit increase in tumour volume, the log of odds of having higher stage increases by 6.3180273
  
  #the slope(i.e. the coefficient) in the regression equation is equal to
  #the logarithm of the odds ratio, so we have to exponentiate the coef to obtain OR
  exp(coef(logistic.model))#odds ration
  exp(cbind(OR = coef(logistic.model), ci))#proportional odds ratio
  #for a 1 unit increase in tumour volume, the odds of moving from one stage to another
  #are 554.4781 times greater, given that the other variables in the model are held constant
  #-the tumour with high volume have 554.4781 times the odds of being very likely to be classified as higher stages
  #compared to students who do not have a lower volume
  exp(-coef(logistic.model))
  #tumours which have a higher tumour volume have a 100-1.8=98.2% lower odds of being less 
  #likely to correspond to a higher stage.
  
  xtumourvol<-seq(min(tumour.vol.rescale),max(tumour.vol.rescale),by=0.01)
  logimd<-polr(stage.MRI.rescalef~tumour.vol.rescale,Hess=TRUE)
  predprob<-predict(logimd,newdata=data.frame(tumour.vol.rescale=xtumourvol))
  #predict.lrm; orm=ordinal logistic model
  #predict.clm
  #type = "response" gives the predicted probabilities
  #type="terms" by default
  png("ggplot_LogitCurve_EC_MRIstage_byTumourVol412.png",width=16,height=8,units='in',res=300)
  plot(tumour.vol.rescale, stage.MRI.rescalef, pch = 16, xlab = "Tumour Volume (mm^3)", ylab = "MRI Stage",main="Logistic Regression Model",col.main="purple")
  lines(xtumourvol, predprob,col="purple")
  dev.off()
  
  nvoxels<-corrected.data$N.Voxels[volpop]
  nvoxels.rescale<-(nvoxels-min(nvoxels))/(max(nvoxels)-min(nvoxels))
  levels(stage.MRI.rescalef)<-levels(stage.MRI)
  
  age.scale<-corrected.data$AgeAtDiagnosis[volpop]
  age.rescale<-(age.scale-min(age.scale))/(max(age.scale)-min(age.scale))
  
  stripchart(tumour.vol.rescale~stage.MRI.rescalef,method="jitter")
  
  boxplot(tumour.vol.rescale~stage.MRI.rescalef)
  #boxplot(quantitative/num.var~qualitative/categ.var)
  
  library(ggplot2)
  wdata<-data.frame(stage.MRI.rescalef,tumour.vol.rescale,nvoxels.rescale,age.rescale)
  png("ggplot_boxplots_EC_MRIstage_byTumourVol412_wnames.png",width=16,height=8,units='in',res=300)
  ggplot(wdata,aes(x = stage.MRI.rescalef, y = tumour.vol.rescale, fill = stage.MRI.rescalef))+   geom_boxplot(size = .75) + labs(x="MRI Stage",y="Tumour Volume (mm^3)",title="Using the Ordinal Logistic Model to predict outcome (from values in the interval where predictor is defined)",fill="MRI Stage")
  dev.off()
  
  library(MASS)
  pr<-profile(logimd)
  plot(pr)
  
  #predict.clm
  xnvoxels<-seq(min(nvoxels.rescale),max(nvoxels.rescale),by=0.01)
  library(brant)
  brant(logimd)
  brant(logistic.model)
  #logimd,newdata=data.frame(tumour.vol.rescale=xtumourvol))
  #library(ordinal)
  predprobs1<-predict(logimd, newdata=data.frame(tumour.vol.rescale))
  png("ggplot_predicted_EC_MRIstage_byTumourVol412.png",width=16,height=8,units='in',res=300)
  ggplot(wdata, aes(x = tumour.vol.rescale, y = predprobs1,colour = stage.MRI.rescalef, group=1)) +geom_line() + labs(x="Tumour Volume (mm^3 but scaled to normal)",y="Probability",title="Model Predictions for all the predictor's values",fill="MRI Stage")#+ facet_grid(nvoxels.rescale~., labeller="label_both",group=1)
  dev.off()
  #Using the Ordinal Logistic Model to predict outcome (from all the predictor's values)"
  #Predictions of the Model taking into consideration all the predictor's values/all the values of the predictor
  # png("ggplot_predicted_EC_MRIstage_byTumourVol412_facettedByAgeGroups.png",width=16,height=8,units='in',res=300)
  # less45<-length(which(corrected.data$AgeAtDiagnosis[volpop]<45))
  # inbetween<-length(which(corrected.data$AgeAtDiagnosis[volpop]>=45&corrected.data$AgeAtDiagnosis[volpop]<=70))
  # more70<-length(which(corrected.data$AgeAtDiagnosis[volpop]>70))
  # age.groups<-c(rep("<45",less45),rep("45-70",inbetween),rep(">70",more70))
  # ggplot(cbind(wdata,age.groups), aes(x = tumour.vol.rescale, y = predprobs, colour = stage.MRI.rescalef,group=1)) +geom_line() + facet_grid(~age.groups)+ labs(x="Tumour Volume (mm^3 but scaled to normal)",y="Probability",title="Predictions of the Ordinal Logistic Model Facetted by Age",fill="MRI Stage")#can only do this faceting with categ variables
  # dev.off()
  
  plot(tumour.vol.rescale,nvoxels.rescale)
  plot(tumour.vol,nvoxels)
  #linear relationship, so we can use pearson's correlation
  cor(tumour.vol,nvoxels)
  cor(tumour.vol.rescale,nvoxels.rescale)
  #perfect positive correlation
  library(ggpubr)
  png("correlation_tumourVol&Nvoxels-corrCoef.png",width=16,height=8,units='in',res=300)
  ggplot(wdata, aes(x=nvoxels, y=tumour.vol)) + 
    geom_point()+
    geom_smooth(method=lm, se=FALSE) + stat_cor(method="pearson")+ labs(x="Number of Voxels",y=expression("Tumour Volume ("~mm^3~")"),title="Relationship between Tumour Volume and Number of Voxels (n = 412)")+ scale_y_continuous(labels = scales::comma,breaks = seq(0, 800000, by = 100000))+ scale_x_continuous(labels = scales::comma,breaks = seq(0, 600000, by = 100000))
  dev.off()
  #aes(label = paste("Pearson correlation R =",..r..,..p.label..,sep = " "))
  
  plot(tumour.vol.rescale,age.rescale)
  plot(tumour.vol,agep)
  cor(tumour.vol.rescale,age.rescale,method="spearman")
  #weak, almost no correlation
  cor(tumour.vol,agep,method="spearman")
  library(ggpubr)
  #pval not significant (neither linear, nor monotonic relationship)
  png("correlation_tumourVol&age.png",width=16,height=8,units='in',res=300)
  ggplot(wdata, aes(x=agep, y=tumour.vol)) + 
    geom_point()+
    geom_smooth(method=lm, se=FALSE)+stat_compare_means(method = "wilcox.test",paired = TRUE)+ labs(x="Age (years)",y="Tumour Volume (mm^3)",title="Relationship between Tumour Volume and Age at Diagnosis")+ scale_y_continuous(labels = scales::comma,breaks = seq(0, 800000, by = 100000))
  dev.off()
  wilcox.test(agep, tumour.vol, paired = TRUE, alternative = "two.sided")
  #we reject the null hypothesis
  #agep and tumour.vol are statistically different
  
  #geom_smooth(method=lm, se=FALSE) +
  #geom_hline(yintercept=mean(tumour.vol))
  #correlation between predictors and outcome
  #already did this for the main predictor 
  #now for the secondary ones for which we adjust the model
  boxplot(age.rescale~stage.MRI.rescalef)
  png("correlation_tumourVol&MRIstage-normalized.png",width=16,height=8,units='in',res=300)
  ggp1<-ggplot(wdata,aes(x = stage.MRI.rescalef, y = tumour.vol.rescale, fill = stage.MRI.rescalef))+   geom_boxplot(size = .75) + labs(x="MRI Stage",y="Tumour Volume (normalized)")+theme(legend.position = "none",plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))#+fill="MRI Stage"
  #,title="Relationship between Tumour Volume and MRI Stage"
  stripchart(age.rescale~stage.MRI.rescalef,method="jitter")
  ggp2<-ggplot(wdata,aes(x = stage.MRI.rescalef, y = tumour.vol.rescale, color = stage.MRI.rescalef))+   geom_jitter(position=position_jitter(0.2)) + labs(x="MRI Stage",y="Tumour Volume (normalized)")+theme(legend.position = "none",plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))#+scale_color_viridis(discrete = TRUE, option = "D")+scale_fill_viridis(discrete = TRUE) +theme_minimal() +theme(legend.position = "bottom")#+theme(legend.position = "none")
  #scale_color_manual(values =  c(viridis(8),"orange"))+
  grid.arrange(ggp1,ggp2, ncol = 2, top = "Relationship between Tumour Volume and MRI Stage\n")
  dev.off()
  
  png("correlation_tumourVol&MRIstage.png",width=16,height=8,units='in',res=300)
  ggp1<-ggplot(wdata,aes(x = stage.MRI, y = tumour.vol, fill = stage.MRI))+   geom_boxplot(size = .75) + scale_y_continuous(labels = scales::comma,breaks = seq(0, 800000, by = 100000))+labs(x="MRI Stage",y="Tumour Volume (mm^3)")+theme(legend.position = "none",plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
  ggp2<-ggplot(wdata,aes(x = stage.MRI, y = tumour.vol, color = stage.MRI))+   geom_jitter(position=position_jitter(0.2))+ scale_y_continuous(labels = scales::comma,breaks = seq(0, 800000, by = 100000)) + labs(x="MRI Stage",y="Tumour Volume (mm^3)")+theme(legend.position = "none",plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
  grid.arrange(ggp1,ggp2, ncol = 2, top = "Relationship between Tumour Volume and MRI Stage\n")
  dev.off()
  
  dfxy<-data.frame(stage.MRI,agep,tumour.vol)
  png("corr_byTumourVol&MRIstage1.png",width=16,height=8,units='in',res=300)
  ggboxplot(dfxy,x="stage.MRI",y="tumour.vol",color="stage.MRI",add="jitter",ggtheme = theme_bw())+ scale_y_continuous(labels = scales::comma,breaks = seq(0, 800000, by = 100000))+labs(x="MRI Stage",y=expression("Tumour Volume ("~mm^3~")"),title="Relationship between Tumour Volume and MRI Stage (n = 412)")+theme(legend.position = "none",text=element_text(family="Arial"))
  dev.off()
  
  is_extreme_outlier <- function(x) {
    #return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
    return(x < quantile(x, 0.25) - 3 * IQR(x) | x > quantile(x, 0.75) + 3 * IQR(x))  
  }
  is_different<-function(x){
    return(x>300000)
      
  }
  for(i in 1:length(tumour.vol))
    if(tumour.vol[i]>300000)
      print(tumour.vol[i])
  #library(ggsignif)
  ll<-list()
  ll["IA"]<-length(which(stage.MRI=="IA"))
  ll["IB"]<-length(which(stage.MRI=="IB"))
  ll["II"]<-length(which(stage.MRI=="II"))
  ll["IIIA"]<-length(which(stage.MRI=="IIIA"))
  ll["IIIB"]<-length(which(stage.MRI=="IIIB"))
  ll["IIIC1"]<-length(which(stage.MRI=="IIIC1"))
  ll["IIIC2"]<-length(which(stage.MRI=="IIIC2"))
  ll["IV"]<-length(which(stage.MRI=="IV"))
  # factor(stage.MRI, labels = paste0("n = ",length(levels(stage.MRI))))
  
  lb<-paste0(levels(stage.MRI),"\n","n = ",ll[levels(stage.MRI)])
  png("corr_byTumourVol&MRIstage.png",width=16,height=8,units='in',res=300)
  ggplot(dfxy,aes(stage.MRI,tumour.vol))+geom_boxplot(outlier.shape = NA)+stat_compare_means(paired=TRUE,label.x=7.2,label.y=700000)+
    geom_jitter(aes(color=stage.MRI),position=position_jitter(0.2))+scale_y_continuous(labels = scales::comma,breaks = seq(0, 800000, by = 100000))+labs(x="MRI Stage",y=expression("Tumour Volume ("~mm^3~")"),title="Relationship between Tumour Volume and MRI Stage (n = 412 in total)")+theme(legend.position = "none")+
    geom_point(shape=ifelse(tumour.vol>300000,8,NA))+#scale_shape_manual(values=tumour.vol,breaks=is_different)
  #+ geom_signif(comparisons = list(c("IA","IB","II","IIIA","IIIB","IIIC1","IIIC2","IV")),annotations = c("*"))
    geom_hline(yintercept=300000)+
    scale_x_discrete(labels=lb)
    
  dev.off()
  
  
  # png("correlation_age&MRIstage.png",width=16,height=8,units='in',res=300)
  # ggp1<-ggplot(wdata,aes(x = stage.MRI, y = agep, fill = stage.MRI))+   geom_boxplot(size = .75) + labs(x="MRI Stage",y="Age (years)")+theme(legend.position = "none",plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
  # ggp2<-ggplot(wdata,aes(x = stage.MRI, y = agep, color = stage.MRI))+   geom_jitter(position=position_jitter(0.2))+ labs(x="MRI Stage",y="Age (years)")+theme(legend.position = "none",plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
  # grid.arrange(ggp1,ggp2, ncol = 2, top = "Relationship between Age at Diagnosis and MRI Stage\n")
  # dev.off()
  
  # png("corr_Byage&MRIstage1.png",width=16,height=8,units='in',res=300)
  # ggplot(wdata,aes(stage.MRI,agep))+geom_boxplot(outlier.shape = NA)+
  #   geom_jitter(aes(color=stage.MRI),position=position_jitter(0.2))+labs(x="MRI Stage",y="Age (years)",title="Relationship between Age at Diagnosis and MRI Stage")+theme(legend.position = "none")
  # #,position=position_dodge(width=0.5))
  # dev.off()
  #improvement:
  png("corr_Byage&MRIstage.png",width=16,height=8,units='in',res=300)
  ggplot(dfxy,aes(stage.MRI,agep))+geom_boxplot(outlier.shape = NA)+stat_compare_means(paired=TRUE,label.x=7.35,label.y=100)+
    geom_jitter(aes(color=stage.MRI),position=position_jitter(0.2))+labs(x="MRI Stage",y="Age (years)",title="Relationship between Age at Diagnosis and MRI Stage (n = 412 in total)")+theme(legend.position = "none")+
    scale_x_discrete(labels=lb)#+scale_y_continuous(breaks=seq(0,100,by=20),limits=c(0,100))
  dev.off()
  # png("corr_Byage&MRIstage2.png",width=16,height=8,units='in',res=300)
  # ggboxplot(dfxy,x="stage.MRI",y="agep",color="stage.MRI",add="jitter",ggtheme = theme_bw())+labs(x="MRI Stage",y="Age (years)",title="Relationship between Age at Diagnosis and MRI Stage (n = 412)")+theme(legend.position = "none")
  # dev.off()
  # png("corr_Byage&MRIstage_dots.png",width=16,height=8,units='in',res=300)
  # ggboxplot(dfxy,x="stage.MRI",y="agep",color="stage.MRI",ggtheme = theme_bw())+labs(x="MRI Stage",y="Age (years)",title="Relationship between Age at Diagnosis and MRI Stage")+theme(legend.position = "none")+
  # stat_compare_means(label="p.signif",method = "wilcox.test")
  # dev.off()
  
  xtabs(~grade.histo.rescalef+stage.MRI.rescalef[-228])
  prop.table(xtabs(~grade.histo.rescalef+stage.MRI.rescalef[-228]),1)
  
  library(MASS)
  library(ordinal)
  clm(stage.MRI~tumour.vol+agep+nvoxels,Hess=TRUE,data=corrected.data)
  #clm function suggests scalling
  olm1<-polr(stage.MRI.rescalef~tumour.vol.rescale+age.rescale,Hess=TRUE,data=corrected.data)
  olm2<-clm(stage.MRI.rescalef~tumour.vol.rescale+age.rescale,Hess=TRUE,data=corrected.data)
  
  summary(logimd)
  summary(olm1)
  anova(logimd,olm1,test="Chisq")
  anova(logimd,olm1)
  olm<-polr(stage.MRI.rescalef~tumour.vol.rescale+age.rescale+nvoxels.rescale,Hess=TRUE,data=corrected.data)
  anova(olm,olm1)#pval not sign, so we stick with olm1
  
  ctable1 <- coef(summary(olm1))
  #table of the regression coefficients
  pval1 <- pnorm(abs(ctable1[, "t value"]), lower.tail = FALSE) * 2
  ctable1 <- cbind(ctable1, "p value" = pval1)
  #advanced stages have the largest absolute values of the coefficients;
  #all have p<0.05, meaning that they are statistically significant
  #we can check this with confindence intervals as well
  ci1 <- confint(olm1)
  confint.default(olm1) # CIs assuming normality
  #the conf.int does not contain 0, so the parameter estimate is stat.sign.
  ##for a unit increase in tumour volume, the log of odds of having higher stage increases by 6.3180273
  #the estimated effect is smaller=>this suggests confounding by age
  
  #the slope(i.e. the coefficient) in the regression equation is equal to
  #the logarithm of the odds ratio, so we have to exponentiate the coef to obtain OR
  exp(coef(olm1))#odds ratio
  exp(cbind(OR = coef(olm1), ci))#proportional odds ratio
  #for a 1 unit increase in tumour volume, the odds of moving from one stage to another
  #are 530 times greater, given that the other variables in the model are held constant
  #-the tumour with high volume have 530 times the odds of being very likely to be classified as higher stages
  #compared to students who do not have a lower volume
  exp(-coef(olm1))
  #tumours which have a higher tumour volume have a 100-1.8=98.2% lower odds of being less 
  #likely to correspond to a higher stage.
  library(brant)
  brant(olm1)
  
  #stage.MRI.rescalef~tumour.vol.rescale+age.rescale
  #nwdata<-data.frame(tumour.vol.rescale=xtumourvol,age.rescale=xage)
  plot(stage.MRI.rescalef~tumour.vol.rescale+age.rescale,data=wdata,col="red")
  #lines(stage.MRI.rescalef~tumour.vol.rescale+age.rescale, data=lnewdat ,col="blue",lwd=2)
  
  library(ordinal)
  xage<-seq(min(age.rescale),max(age.rescale),by=0.01)
  predsprob<-predict(olm1,newdata=data.frame(tumour.vol.rescale=xtumourvol,age.rescale=xage))
  png("predictedMRIstageCateg_byTumourVol&Age412.png",width=15,height=8,units='in',res=300)
  plot(tumour.vol.rescale, stage.MRI.rescalef, pch = 16, xlab = "Tumour Volume (normalized)", ylab = "MRI Stage (normalized)",main="The MRI Stage Predicted by Two Ordinal Regression Models (n = 412)")
  #main="Logistic Regression Model Adjusted by Age",col.main="purple"
  lines(xtumourvol, predprob,col="purple")
  lines(xtumourvol,predsprob,col="blue",lty="dashed")
  #par(xpd=TRUE)
  text(0.45,6,"p = 0.007187584 for chi-squared test",pos=4)
  legend("bottomright",title="Models together with the p-value for the Lipsitz test",legend=c("with tumour volume alone,\np = 3.84e-09\n", "with tumour volume adjusted for age,\np = 7.034e-05\n"),
         col=c("purple", "blue"), lty=1:2, cex=1)
  axis(2, at=1:8, labels=levels(stage.MRI),pos=-0.055,tick=FALSE)
  dev.off()
  
  library(generalhoslem)
  lt1<-lipsitz.test(olm1)
  lt1$p.value
  lt2<-lipsitz.test(logimd)
  lt2$p.value
  
  library(ggeffects)
  plot(ggpredict(olm1,"tumour.vol.rescale"))
  plot(ggpredict(olm1,"age.rescale"))
  
  logit.df<-data.frame(xagec=seq.int(min(agep),max(agep),len=412),xtumvolc=seq(min(tumour.vol), max(tumour.vol),len=412))
  # library(rms)
  # library(glm)
  # library(polr)
  # library(MASS)
  library(lrm)
  predictedp<-predict(olm1,logit.df,type="probs")
  plot(predictedp,data=corrected.data,pch = 16, xlab = "Tumour Volume (normalized)", ylab = "MRI Stage (normalized)",main="The Probabilities Predicted by Two Ordinal Regression Models")#Logistic Regression Curves
  #provided by the function predict
  #lines(~, logit.df, col="red", lwd=2)
  
  require(ggiraph)
  require(ggiraphExtra)
  require(plyr)
  ggPredict(olm1)
  
  
  #error#effect_plot(olm2, pred = tumour.vol.rescale)
  #pred should be part of the model#effect_plot(olm1, pred = tumour.vol)
  sjPlot::plot_model(olm1)
  
  prob<-predict(logimd, newdata=data.frame(tumour.vol.rescale),type="probs")
  png("ggplot_predictedProb_EC_MRIstageIA_byTumourVol412.png",width=16,height=8,units='in',res=300)
  ggplot(wdata, aes(x = tumour.vol.rescale, y = prob[,1])) +geom_line() + labs(x="Tumour Volume (normalized)",y="Probability",title="Probabilities (that tumour vol will fall within stage IA) predicted by the Ordinal Logistic Model",fill="MRI Stage")#can only do this faceting with categ variables
  dev.off()
  #+ facet_grid(~age.groups)+
  png("ggplot_predictedProb_EC_MRIstageIB_byTumourVol412.png",width=16,height=8,units='in',res=300)
  ggplot(wdata, aes(x = tumour.vol.rescale, y = prob[,2])) +geom_line() + labs(x="Tumour Volume (normalized)",y="Probability",title="Probabilities (that tumour vol will fall within stage IB) predicted by the Ordinal Logistic Model",fill="MRI Stage")#can only do this faceting with categ variables
  dev.off()
  png("ggplot_predictedProb_EC_MRIstageII_byTumourVol412.png",width=16,height=8,units='in',res=300)
  ggplot(wdata, aes(x = tumour.vol.rescale, y = prob[,3])) +geom_line() + labs(x="Tumour Volume (normalized)",y="Probability",title="Probabilities (that tumour vol will fall within stage II) predicted by the Ordinal Logistic Model",fill="MRI Stage")#can only do this faceting with categ variables
  dev.off()
  png("ggplot_predictedProb_EC_MRIstageIIIA_byTumourVol412.png",width=16,height=8,units='in',res=300)
  ggplot(wdata, aes(x = tumour.vol.rescale, y = prob[,4])) +geom_line() + labs(x="Tumour Volume (normalized)",y="Probability",title="Probabilities (that tumour vol will fall within stage IIIA) predicted by the Ordinal Logistic Model",fill="MRI Stage")#can only do this faceting with categ variables
  dev.off()
  png("ggplot_predictedProb_EC_MRIstageIIIB_byTumourVol412.png",width=16,height=8,units='in',res=300)
  ggplot(wdata, aes(x = tumour.vol.rescale, y = prob[,5])) +geom_line() + labs(x="Tumour Volume (normalized)",y="Probability",title="Probabilities (that tumour vol will fall within stage IIIB) predicted by the Ordinal Logistic Model",fill="MRI Stage")#can only do this faceting with categ variables
  dev.off()
  png("ggplot_predictedProb_EC_MRIstageIIIC1_byTumourVol412.png",width=16,height=8,units='in',res=300)
  ggplot(wdata, aes(x = tumour.vol.rescale, y = prob[,6])) +geom_line() + labs(x="Tumour Volume (normalized)",y="Probability",title="Probabilities (that tumour vol will fall within stage IIIC1) predicted by the Ordinal Logistic Model",fill="MRI Stage")#can only do this faceting with categ variables
  dev.off()
  png("ggplot_predictedProb_EC_MRIstageIIIC2_byTumourVol412.png",width=16,height=8,units='in',res=300)
  ggplot(wdata, aes(x = tumour.vol.rescale, y = prob[,7])) +geom_line() + labs(x="Tumour Volume (normalized)",y="Probability",title="Probabilities (that tumour vol will fall within stage IIIC2) predicted by the Ordinal Logistic Model",fill="MRI Stage")#can only do this faceting with categ variables
  dev.off()
  png("ggplot_predictedProb_EC_MRIstageIV_byTumourVol412.png",width=16,height=8,units='in',res=300)
  ggplot(wdata, aes(x = tumour.vol.rescale, y = prob[,8])) +geom_line() + labs(x="Tumour Volume (normalized)",y="Probability",title="Probabilities (that tumour vol will fall within stage IV) predicted by the Ordinal Logistic Model",fill="MRI Stage")#can only do this faceting with categ variables
  dev.off()
  
  png("ggplot_predictedProb_EC_MRIstageIA_bymm3TumourVol412.png",width=16,height=8,units='in',res=300)
  ggplot(wdata, aes(x = tumour.vol, y = prob[,1])) +geom_line() + labs(x="Tumour Volume (mm^3)",y="Probability",title="Probabilities (that tumour vol will fall within stage IA) predicted by the Ordinal Logistic Model",fill="MRI Stage")#can only do this faceting with categ variables
  dev.off()
  #+ facet_grid(~age.groups)+
  png("ggplot_predictedProb_EC_MRIstageIB_bymm3TumourVol412.png",width=16,height=8,units='in',res=300)
  ggplot(wdata, aes(x = tumour.vol, y = prob[,2])) +geom_line() + labs(x="Tumour Volume (mm^3)",y="Probability",title="Probabilities (that tumour vol will fall within stage IB) predicted by the Ordinal Logistic Model",fill="MRI Stage")#can only do this faceting with categ variables
  dev.off()
  png("ggplot_predictedProb_EC_MRIstageII_bymm3TumourVol412.png",width=16,height=8,units='in',res=300)
  ggplot(wdata, aes(x = tumour.vol, y = prob[,3])) +geom_line() + labs(x="Tumour Volume (mm^3)",y="Probability",title="Probabilities (that tumour vol will fall within stage II) predicted by the Ordinal Logistic Model",fill="MRI Stage")#can only do this faceting with categ variables
  dev.off()
  png("ggplot_predictedProb_EC_MRIstageIIIA_bymm3TumourVol412.png",width=16,height=8,units='in',res=300)
  ggplot(wdata, aes(x = tumour.vol, y = prob[,4])) +geom_line() + labs(x="Tumour Volume (mm^3)",y="Probability",title="Probabilities (that tumour vol will fall within stage IIIA) predicted by the Ordinal Logistic Model",fill="MRI Stage")#can only do this faceting with categ variables
  dev.off()
  png("ggplot_predictedProb_EC_MRIstageIIIB_bymm3TumourVol412.png",width=16,height=8,units='in',res=300)
  ggplot(wdata, aes(x = tumour.vol, y = prob[,5])) +geom_line() + labs(x="Tumour Volume (mm^3)",y="Probability",title="Probabilities (that tumour vol will fall within stage IIIB) predicted by the Ordinal Logistic Model",fill="MRI Stage")#can only do this faceting with categ variables
  dev.off()
  png("ggplot_predictedProb_EC_MRIstageIIIC1_bymm3TumourVol412.png",width=16,height=8,units='in',res=300)
  ggplot(wdata, aes(x = tumour.vol, y = prob[,6])) +geom_line() + labs(x="Tumour Volume (mm^3)",y="Probability",title="Probabilities (that tumour vol will fall within stage IIIC1) predicted by the Ordinal Logistic Model",fill="MRI Stage")#can only do this faceting with categ variables
  dev.off()
  png("ggplot_predictedProb_EC_MRIstageIIIC2_bymm3TumourVol412.png",width=16,height=8,units='in',res=300)
  ggplot(wdata, aes(x = tumour.vol, y = prob[,7])) +geom_line() + labs(x="Tumour Volume (mm^3)",y="Probability",title="Probabilities (that tumour vol will fall within stage IIIC2) predicted by the Ordinal Logistic Model",fill="MRI Stage")#can only do this faceting with categ variables
  dev.off()
  png("ggplot_predictedProb_EC_MRIstageIV_bymm3TumourVol412.png",width=16,height=8,units='in',res=300)
  ggplot(wdata, aes(x = tumour.vol, y = prob[,8])) +geom_line() + labs(x="Tumour Volume (mm^3)",y="Probability",title="Probabilities (that tumour vol will fall within stage IV) predicted by the Ordinal Logistic Model",fill="MRI Stage")#can only do this faceting with categ variables
  dev.off()
  
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(ggtext)
  library(viridis)
  png("predictedProbVaryingByMRIstage.png",width=16,height=16,units='in',res=300)
  #t1 <- textGrob(expression("Concentration of " * phantom(bold("affluence")) * "and" * phantom(bold("poverty")) * " nationwide"),
                 #gp = gpar(col = "black"))
  
  #t2 <- textGrob(expression(phantom("Concentration of ") * bold("affluence") * phantom(" and poverty nationwide")),
                 #gp = gpar(col = "#EEB422"))
  
  #t3 <- textGrob(expression(phantom("Concentration of affluence and ") * bold("poverty") * phantom(" nationwide")),
                #gp = gpar(col = "#238E68"))
  t1 <- textGrob(expression(bold("Probability that tumour volume falls within ") * phantom("stage IA")),x=0.475,y=1, gp = gpar(col = "black"))
  t2 <- textGrob(expression(phantom("Probability that tumour volume falls within ") * bold("stage IA")),x=0.5,y=1, gp = gpar(col = viridis(8)[1]))
  p1<-ggplot(wdata, aes(x = tumour.vol, y = prob[,1])) +geom_line()+  annotation_custom(grobTree(t1,t2))+ scale_x_continuous(labels = scales::comma,breaks = seq(0, 800000, by = 100000)) +coord_cartesian(clip = "off") +labs(x=expression("Tumour Volume ("~mm^3~")"),y="Probability",fill="MRI Stage") +theme_minimal()+theme(plot.background = element_rect(color = viridis(8)[1], size = 1),legend.position = 'none',plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))#+annotation_custom(grobTree(t1,t2))
  #plot.title=element_text(margin = margin(b = -0.5)),                                                                                                                                                                                                                                                      
  #grid.rect(width = .98, height = .98, gp = gpar(lwd = 2, col = "blue", fill = NA))
  #+scale_y_continuous(breaks=seq(0,1,by=0.1))#limits = c(0, 800000)
  #annotate(geom = "point", x = min(tumour.vol), y = 0, size = 1, color = 'red')+annotate(geom = "point", x = max(tumour.vol), y = 0, size = 1, color = 'red')+
  #gp = gpar(col = "#238E68"))
  t1 <- textGrob(expression(bold("Probability that tumour volume falls within ") * phantom("stage IB")),x=0.475,y=1, gp = gpar(col = "black"))
  t2 <- textGrob(expression(phantom("Probability that tumour volume falls within ") * bold("stage IB")),x=0.5,y=1, gp = gpar(col = viridis(8)[2]))
  p2<-ggplot(wdata, aes(x = tumour.vol, y = prob[,2])) +geom_line()+  annotation_custom(grobTree(t1,t2))+ scale_x_continuous(labels = scales::comma,breaks = seq(0, 800000, by = 100000)) +coord_cartesian(clip = "off") +labs(x=expression("Tumour Volume ("~mm^3~")"),y="Probability",fill="MRI Stage") +theme_minimal()+theme(plot.background = element_rect(color = viridis(8)[2], size = 1),legend.position = 'none',plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
  t1 <- textGrob(expression(bold("Probability that tumour volume falls within ") * phantom("stage II")),x=0.475,y=1, gp = gpar(col = "black"))
  t2 <- textGrob(expression(phantom("Probability that tumour volume falls within ") * bold("stage II")),x=0.5,y=1, gp = gpar(col = viridis(8)[3]))
  p3<-ggplot(wdata, aes(x = tumour.vol, y = prob[,3])) +geom_line()+  annotation_custom(grobTree(t1,t2))+ scale_x_continuous(labels = scales::comma,breaks = seq(0, 800000, by = 100000)) +coord_cartesian(clip = "off") +labs(x=expression("Tumour Volume ("~mm^3~")"),y="Probability",fill="MRI Stage") +theme_minimal()+theme(plot.background = element_rect(color = viridis(8)[3], size = 1),legend.position = 'none',plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
  t1 <- textGrob(expression(bold("Probability that tumour volume falls within ") * phantom("stage IIIA")),x=0.475,y=1, gp = gpar(col = "black"))
  t2 <- textGrob(expression(phantom("Probability that tumour volume falls within ") * bold("stage IIIA")),x=0.5,y=1, gp = gpar(col = viridis(8)[4]))
  p4<-ggplot(wdata, aes(x = tumour.vol, y = prob[,4])) +geom_line()+  annotation_custom(grobTree(t1,t2))+ scale_x_continuous(labels = scales::comma,breaks = seq(0, 800000, by = 100000)) +coord_cartesian(clip = "off") +labs(x=expression("Tumour Volume ("~mm^3~")"),y="Probability",fill="MRI Stage") +theme_minimal()+theme(plot.background = element_rect(color = viridis(8)[4], size = 1),legend.position = 'none',plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
  t1 <- textGrob(expression(bold("Probability that tumour volume falls within ") * phantom("stage IIIB")),x=0.475,y=1, gp = gpar(col = "black"))
  t2 <- textGrob(expression(phantom("Probability that tumour volume falls within ") * bold("stage IIIB")),x=0.5,y=1, gp = gpar(col = viridis(8)[5]))
  p5<-ggplot(wdata, aes(x = tumour.vol, y = prob[,5])) +geom_line()+  annotation_custom(grobTree(t1,t2))+ scale_x_continuous(labels = scales::comma,breaks = seq(0, 800000, by = 100000)) +coord_cartesian(clip = "off") +labs(x=expression("Tumour Volume ("~mm^3~")"),y="Probability",fill="MRI Stage") +theme_minimal()+theme(plot.background = element_rect(color = viridis(8)[5], size = 1),legend.position = 'none',plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
  t1 <- textGrob(expression(bold("Probability that tumour volume falls within ") * phantom("stage IIIC1")),x=0.475,y=1, gp = gpar(col = "black"))
  t2 <- textGrob(expression(phantom("Probability that tumour volume falls within ") * bold("stage IIIC1")),x=0.5,y=1, gp = gpar(col = viridis(8)[6]))
  p6<-ggplot(wdata, aes(x = tumour.vol, y = prob[,6])) +geom_line()+  annotation_custom(grobTree(t1,t2))+ scale_x_continuous(labels = scales::comma,breaks = seq(0, 800000, by = 100000)) +coord_cartesian(clip = "off") +labs(x=expression("Tumour Volume ("~mm^3~")"),y="Probability",fill="MRI Stage") +theme_minimal()+theme(plot.background = element_rect(color = viridis(8)[6], size = 1),legend.position = 'none',plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
  t1 <- textGrob(expression(bold("Probability that tumour volume falls within ") * phantom("stage IIIC2")),x=0.475,y=1, gp = gpar(col = "black"))
  t2 <- textGrob(expression(phantom("Probability that tumour volume falls within ") * bold("stage IIIC2")),x=0.5,y=1, gp = gpar(col = viridis(8)[7]))
  p7<-ggplot(wdata, aes(x = tumour.vol, y = prob[,7])) +geom_line()+  annotation_custom(grobTree(t1,t2))+ scale_x_continuous(labels = scales::comma,breaks = seq(0, 800000, by = 100000)) +coord_cartesian(clip = "off") +labs(x=expression("Tumour Volume ("~mm^3~")"),y="Probability",fill="MRI Stage") +theme_minimal()+theme(plot.background = element_rect(color = viridis(8)[7], size = 1),legend.position = 'none',plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
  t1 <- textGrob(expression(bold("Probability that tumour volume falls within ") * phantom("stage IV")),x=0.475,y=1, gp = gpar(col = "black"))
  t2 <- textGrob(expression(phantom("Probability that tumour volume falls within ") * bold("stage IV")),x=0.5,y=1, gp = gpar(col = "orange"))
  p8<-ggplot(wdata, aes(x = tumour.vol, y = prob[,8])) +geom_line()+  annotation_custom(grobTree(t1,t2))+ scale_x_continuous(labels = scales::comma,breaks = seq(0, 800000, by = 100000)) +coord_cartesian(clip = "off") +labs(x=expression("Tumour Volume ("~mm^3~")"),y="Probability",fill="MRI Stage") +theme_minimal()+theme(plot.background = element_rect(color = "orange", size = 1),legend.position = 'none',plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
  #p2<-ggplot(wdata, aes(x = tumour.vol, y = prob[,2])) +geom_line() + labs(x="Tumour Volume (mm^3)",y="Probability",title="Probability that tumour volume falls within stage IB",fill="MRI Stage") + scale_x_continuous(labels = scales::comma,breaks = seq(0, 800000, by = 100000))
  #p3<-ggplot(wdata, aes(x = tumour.vol, y = prob[,3])) +geom_line() + labs(x="Tumour Volume (mm^3)",y="Probability",title="Probability that tumour volume falls within stage II",fill="MRI Stage") + scale_x_continuous(labels = scales::comma,breaks = seq(0, 800000, by = 100000))
  #p4<-ggplot(wdata, aes(x = tumour.vol, y = prob[,4])) +geom_line() + labs(x="Tumour Volume (mm^3)",y="Probability",title="Probability that tumour volume falls within stage IIIA",fill="MRI Stage") + scale_x_continuous(labels = scales::comma,breaks = seq(0, 800000, by = 100000))
  #p5<-ggplot(wdata, aes(x = tumour.vol, y = prob[,5])) +geom_line() + labs(x="Tumour Volume (mm^3)",y="Probability",title="Probability that tumour volume falls within stage IIIB",fill="MRI Stage") + scale_x_continuous(labels = scales::comma,breaks = seq(0, 800000, by = 100000))
  #p6<-ggplot(wdata, aes(x = tumour.vol, y = prob[,6])) +geom_line() + labs(x="Tumour Volume (mm^3)",y="Probability",title="Probability that tumour volume falls within stage IIIC1",fill="MRI Stage") + scale_x_continuous(labels = scales::comma,breaks = seq(0, 800000, by = 100000))
  #p7<-ggplot(wdata, aes(x = tumour.vol, y = prob[,7])) +geom_line() + labs(x="Tumour Volume (mm^3)",y="Probability",title="Probability that tumour volume falls within stage IIIC2",fill="MRI Stage") + scale_x_continuous(labels = scales::comma,breaks = seq(0, 800000, by = 100000))
  #p8<-ggplot(wdata, aes(x = tumour.vol, y = prob[,8])) +geom_line() + labs(x="Tumour Volume (mm^3)",y="Probability",title="Probability that tumour volume falls within stage IV",fill="MRI Stage") + scale_x_continuous(labels = scales::comma,breaks = seq(0, 800000, by = 100000))
  grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, ncol = 2, top = "Probabilities Predicted by the Ordinal Logistic Model with Tumour Volume Alone\n", bottom="n = 412")
  #(that tumour volumes will fall within certain stage)
  dev.off()
  
  library(ordinal)
  predsprobs<-predict(olm1, newdata=data.frame(tumour.vol.rescale,age.rescale),type="probs")
  # png("ggplot_predicted_EC_MRIstage_byTumourVol&Age412_facettedByAgeGroups.png",width=16,height=8,units='in',res=300)
  # ggplot(cbind(wdata,age.groups), aes(x = tumour.vol.rescale, y = predsprobs, colour = stage.MRI.rescalef,group=1)) +geom_line() + facet_grid(~age.groups)+ labs(x="Tumour Volume (normalized)",y="Probability",title="Predictions of the Ordinal Logistic Model Facetted by Age",fill="MRI Stage")#can only do this faceting with categ variables
  # dev.off()
  newdt<-data.frame(tumour.vol,agep)
  newdat <- cbind(newdt, predict(olm1, newdt, type = "probs"))
  library(reshape2)
  lnewdat <- melt(newdat, id.vars = c("agep", "tumour.vol"),
                  variable.name = "stage.MRI", value.name="Probability")
  head(lnewdat)
  age.groups<-vector()
  age.groups[which(agep<45)]<-'<45'
  age.groups[which(agep>=45&agep<=70)]<-'45-70'
  age.groups[which(agep>70)]<-'>70'
  age.group<-age.groups
  
  png("predictedProbsVaryingByMRIstage-facettedByAgeGroups.png",width=16,height=8,units='in',res=300)
  ggplot(cbind(lnewdat,age.group), aes(x = tumour.vol, y = Probability, colour = stage.MRI)) +
    geom_line() + facet_grid(age.group~.,labeller="label_both")+scale_x_continuous(labels = scales::comma,breaks = seq(0, 800000, by = 100000))+labs(x=expression("Tumour Volume ("~mm^3~")"),title="Predictive Probabilities of the Ordinal Logistic Model with Age + Tumour Volume, Facetted by Age Groups (n = 412)",colour="MRI Stage")
  dev.off()
  
  ggplot(cbind(lnewdat,age.group,nvoxels), aes(x = tumour.vol, y = Probability, colour = stage.MRI)) +
    geom_line() + facet_grid(age.group~.,labeller="label_both")+scale_x_continuous(labels = scales::comma,breaks = seq(0, 800000, by = 100000))+labs(x="Tumour Volume (mm^3)",colour="MRI Stage")
  
  library(ordinal)
  coeff<-coef(olm1)
  #rs<-predict(olm1,type="cum.prob")
  fun<-function(x1,x2) coeff[1]*x1+coeff[2]*x2+1.087837
  x1=tumour.vol.rescale
  x2=age.rescale
  y<-fun(x1,x2)
  tbl<-cbind(log(y),predsprobs[,1])
  head(tbl)
  
  logit<-1.087837-coeff[1]*x1-coeff[2]*x2
  odds <- exp(logit)
  probf <- odds / (1 + odds)
  tblf<-cbind(probf,predsprobs[,1])
  head(tblf)
  tail(tblf)
  
  library(rgl)
  #open3d()
  plot3d(logit, col = colorRampPalette(c("blue", "white", "red")), 
         xlab = "Tumour Volume", ylab = "Age", zlab = "Logit",
         xlim = c(200,80000), ylim = c(20, 100),
         aspect = c(1, 1, 0.5))
  
  #y<-mapply(fun,x1,x2)
  #matplot(x1,y=,main="Logistic Regression Equation for stage IA",type="l",lwd=2)
  #ggplot(cbind(wdata,age.groups), aes(x = tumour.vol.rescale, y = predprobs1, colour = stage.MRI.rescalef,group=1)) +geom_line() + facet_grid(~age.groups)+ labs(x="Tumour Volume (mm^3 but scaled to normal)",y="Probability",title="Predictions of the Ordinal Logistic Model Facetted by Age",fill="MRI Stage")
  
  #if we add a var that is very correlated with one of the other predictors,
  #we find that those covariates are no longer stat. sign. for the model
  # olm<-polr(stage.MRI.rescalef~tumour.vol.rescale+age.rescale+nvoxels.rescale,Hess=TRUE,data=corrected.data)
  # clm(stage.MRI.rescalef~tumour.vol.rescale+age.rescale+nvoxels.rescale,Hess=TRUE,data=corrected.data)
  # 
  # ctabl <- coef(summary(olm))
  # #table of the regression coefficients
  # pval <- pnorm(abs(ctabl[, "t value"]), lower.tail = FALSE) * 2
  # ctabl <- cbind(ctabl, "p value" = pval)
  # #advanced stages have the largest absolute values of the coefficients;
  # #all have p<0.05, meaning that they are statistically significant
  # #we can check this with confindence intervals as well
  # CI <- confint(olm)
  # confint.default(olm) # CIs assuming normality
  # #the conf.int does not contain 0, so the parameter estimate is stat.sign.
  
  #facet_wrap(~nvoxels.rescale)
  #facet_grid(x ~ y) will display x*y plots even if some plots are empty
  #facet_wrap(x ~ y) on the other hand, displays only the plots having actual values.
  
  #gg<-ggplot(wdata,aes(x = stage.MRI.rescalef, y = tumour.vol.rescale, fill = stage.MRI.rescalef))+   geom_boxplot(size = .75)+   facet_grid(nvoxels.rescale ~ age.rescale, margins = FALSE) +   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  #different ways of printing
  #print(gg)
  #library(grid)
  #grid.draw(ggplotGrob(gg))
  
  # png("ggplot_boxplots_EC_MRIstage_byTumourVol412_adjusted.png",width=16,height=8,units='in',res=300)
  # ggplot(wdata,aes(x = stage.MRI.rescalef, y = tumour.vol.rescale, fill = stage.MRI.rescalef))+   geom_boxplot(size = .75)+   facet_grid(nvoxels.rescale ~ age.rescale, margins = FALSE) 
  # +   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  # dev.off()
  
  
  #predict type can be "class" or "probs"
  ##
  #pnorm(2.75,lower.tail = FALSE)*2
  #to see if population mean is more/less than a certain value;
  #compute to the right/left
  #the same principles apply to pt function, just that I have to specify degrees of freedom there
  
  #2*pt(q=2.75, df=5, lower.tail=FALSE)#directly from the table
  ##
  
  p.pt<-2*pt(q=ctable[, "t value"], df=length(stage.MRI.rescalef)-1, lower.tail=FALSE)
  cbind(ctable, "p value from t table" = p.pt)
}
