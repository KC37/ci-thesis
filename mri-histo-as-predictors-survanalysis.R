diagnostic.test<-function(){
  MRI.depth<-corrected.data$Depth.of.Myometrial.Invasion[allstudypop]
  histo.depth<-corrected.data$Sup.Deep.2[allstudypop]
  corrected.data$Sup.Deep.2[which(corrected.data$Sup.Deep.2 %in% c("no","none ","endometrium intact","nil confined to endometrium"))]<-"none"
  corrected.data$Sup.Deep.2[grep("sup",corrected.data$Sup.Deep.2,ignore.case=TRUE)]<-"superficial"
  corrected.data$Sup.Deep.2[grep("deep",corrected.data$Sup.Deep.2,ignore.case=TRUE)]<-"deep"
  corrected.data$Depth.of.Myometrial.Invasion[which(grepl("sup",corrected.data$Depth.of.Myometrial.Invasion,ignore.case=TRUE)&!grepl("deep",corrected.data$Depth.of.Myometrial.Invasion))]<-"superficial"
  corrected.data$Depth.of.Myometrial.Invasion[which(grepl("deep",corrected.data$Depth.of.Myometrial.Invasion,ignore.case=TRUE)&!grepl("no",corrected.data$Depth.of.Myometrial.Invasion)&!grepl("superficial",corrected.data$Depth.of.Myometrial.Invasion))]<-"deep"
  #filled.ind<-intersect(allstudypop,which(corrected.data$Sup.Deep.2 %nin% c("","y","no comment") & corrected.data$Depth.of.Myometrial.Invasion %nin% c("no deep invasion","no suspicious mass post TCRE")))
  #minv.histo<-corrected.data$Sup.Deep.2[filled.ind]
  #minv.MRI<-corrected.data$Depth.of.Myometrial.Invasion[filled.ind]
  avail.histo<-intersect(early.st.pop,which(corrected.data$STAGE %in% c("IA","IB")))
  avail.stageI<-intersect(avail.histo,which(corrected.data$Sup.Deep.2 %in% c("deep","none","superficial")))
  histo.notdeep<-intersect(avail.stageI,which(corrected.data$Sup.Deep.2 %in% c("none","superficial","no deep invasion")))
  histo.deep<-intersect(avail.stageI,which(corrected.data$Sup.Deep.2 == "deep"))
  MRI.notdeep<-intersect(avail.stageI,which(corrected.data$Depth.of.Myometrial.Invasion %in% c("superficial","none","no deep invasion")))
  MRI.deep<-intersect(avail.stageI,which(corrected.data$Depth.of.Myometrial.Invasion == "deep"))
  minv.histo<-vector()
  minv.histo[histo.notdeep]<-"no"
  minv.histo[histo.deep]<-"yes"
  minv.MRI<-vector()
  minv.MRI[MRI.notdeep]<-"no"
  minv.MRI[MRI.deep]<-"yes"
  #paired data because for the same patient we assessed depth of invasion based on MRI and histology
  #we perform McNemar chi-squared test for our binary data (deep/no deep invasion)
  #we build the contingency table
  
  minv.data<-table(minv.histo,minv.MRI)
  mcnemar.test(minv.data)
  mcnemar.test(minv.data,correct=FALSE)
  #p-val>0.5>0.05=>we cannot reject the null hypothesis that 
  #the depth of myometrial invasion is the same on MRI and histology 
  
}

MRIagainstHisto<-function()
{
  final.df<-data.frame(ec.surv1[,1],ec.surv1[,2],stage.MRI1,age1,tumourv.groups2,stage.histo1)
  sbs<-which(early.st.pop %in% avail.stageI)
  model11shisto<-coxph(formula = ec.surv1 ~ stage.MRI1 + stage.histo1 + age1 + 
                        tumourv.groups2,subset=sbs)
  model11shistoA<-coxph(ec.surv1 ~ stage.histo1 + age1 + tumourv.groups2,subset=sbs)
  model11shistoB<-coxph(formula = ec.surv1 ~ stage.MRI1 + age1 + tumourv.groups2,
                        subset=sbs)
  png("MRI&Histo_stage-OS_subset.png",width=10,height=8,units='in',res=300)
  ggls<-list(
    "with MRI stage + histological stage + age + tumour.groupingII"=survfit(model11shisto,data=final.df,subset=sbs),
    "with histological stage + age + tumour.groupingII"=survfit(model11shistoA,data=final.df,subset=sbs),
    "with MRI stage + age + tumour.groupingII"=survfit(model11shistoB,data=final.df,subset=sbs))
  ggsc<-ggsurvplot_combine(ggls,censor=FALSE,xlab="OS Time (days)",title="Cox (PH) Regression Model Comparisons by Predicted Survival Curves",legend.title = "Model by predictors", final.df[sbs,])#, palette = "Dark2")
  ggsc$plot<-ggsc$plot+scale_colour_manual(values = c(rgb(1,0,0,1),rgb(0,0,1,0.5),rgb(0.9,0.65,0)))
  grid.arrange(ggsc$plot,bottom="\nn = 263, number of events = 33 (subset: MRI stage = I & histological stage = I)")
  dev.off()
  
  png("MRIvsHistostage-predictorstrength-CoxModelsComparisons_subset.png",width=16,height=8,units='in',res=300)
  phist<-plot_summs(model11shisto,model11shistoA,model11shistoB,colors=c("red","blue","orange"),coefs = c("MRI stage = IB relative to IA"="stage.MRI1IB",
                                                                                                   "Histological Stage = IB relative to IA"="stage.histo1IB"
  ),legend.title="Model by predictors",model.names = c("with MRI stage + histological stage + age + tumour.groupingII","with histological stage + age + tumour.groupingII","with MRI stage + age + tumour.groupingII"))
  phist<-phist+labs(title="Comparing the Primary Predictors Strength (MRI vs Histological Stage):\nCox (PH) Regression Model Comparisons by Beta Coefficients",x="Estimate (Regression Coefficient)",y="Predictor")+theme(text=element_text(family="Arial"))
  grid.arrange(phist,bottom="n = 263, number of events = 33 (subset: MRI stage = I & histological stage = I)")
  dev.off()
  
  table(stage.MRI1[histo1])
  table(stage.histo1[histo1])
  
  diagnostic.test()
}
