hadRadiot<-intersect(radiot.pop,which(radiotherapy=="yes"))
hadNRadiot<-intersect(radiot.pop,which(radiotherapy=="no"))
table(corrected.data$Risk.score[hadRadiot])
table(corrected.data$Risk.score[hadNRadiot])

table(radiotherapy[which(radiotherapy %in% c("yes","no"))])
table(radiotherapy[which(radiotherapy %nin% c("yes","no"))])
table(chemotherapy[which(chemotherapy %in% c("yes","no"))])
table(chemotherapy[which(chemotherapy %nin% c("yes","no"))])
wtreatment<-intersect(which(radiotherapy %in% c("yes","no")),which(chemotherapy %in% c("yes","no")))
radiot.pop.allst<-intersect(wtreatment,allstudypop)
radiot.pop<-intersect(radiot.pop.allst,which(corrected.data$X.2 %in% c("IA","IB")))
radiot.surv<-generate.Surv.obj(radiot.pop)
age<-corrected.data$AgeAtDiagnosis[radiot.pop]
radiot<-radiotherapy[radiot.pop]
chemot<-chemotherapy[radiot.pop]
gradeh<-corrected.data$Grade[radiot.pop]
typeh<-corrected.data$Histology.Type[radiot.pop]
stageh<-corrected.data$STAGE[radiot.pop]
stageRMN<-corrected.data$X.2[radiot.pop]
CRS<-corrected.data$Risk.score[radiot.pop]
summary(coxph(radiot.surv~radiot+chemot+age+gradeh+typeh+stageh+stageRMN))
md1<-coxph(radiot.surv~radiot)
summary(md1)
md2<-coxph(radiot.surv~radiot+chemot)
summary(md2)
anova(md1,md2,test="LRT")
AIC(md1,md2)
BIC(md1,md2)
#md2 remains
md3<-coxph(radiot.surv~radiot+chemot+stageRMN)
summary(md3)
anova(md2,md3,test="LRT")
#md2 remains
md4<-coxph(radiot.surv~radiot+chemot+stageRMN+gradeh)
summary(md4)
anova(md2,md4,test="LRT")
AIC(md2,md4)
BIC(md2,md4)
#md4 remains
md5<-coxph(radiot.surv~radiot+chemot+stageRMN+gradeh+age)
summary(md5)
anova(md4,md5,test="LRT")
AIC(md4,md5)
BIC(md4,md5)
#md5 remains
md6<-coxph(radiot.surv~radiot+chemot+stageRMN+gradeh+age+stageh)
summary(md6)
anova(md5,md6,test="LRT")
#md5 remains

md7<-coxph(radiot.surv~radiot+chemot+CRS)
summary(md7)
anova(md2,md7,test="LRT")
AIC(md2,md7)
BIC(md2,md7)
#md7 remains
md8<-coxph(radiot.surv~radiot+chemot+CRS+age)
summary(md8)
anova(md7,md8,test="LRT")
AIC(md7,md8)
BIC(md7,md8)
#md8 remains

stime<-radiot.surv[,1]
sevent<-radiot.surv[,2]
dat<-data.frame(stime,sevent,radiot.surv,radiot,chemot,stageRMN,gradeh,age,CRS)
#corrected.data[which(corrected.data$PatientID%in%radiot.pop)]

library(survival)
library(ggplot2)
library(survminer)

#p.adjust for FDR adjustment
mainfit<-survfit(radiot.surv~radiot+stageRMN)
gg2<-ggsurvplot(mainfit,dat,risk.table = TRUE,#pval = TRUE,
           title="The Kaplan-Meier Estimate with Radiotherapy + MRI Stage",
           legend.labs=c("radiotherapy = no, MRI stage = IA",
           "radiotherapy = no, MRI stage = IB",
           "radiotherapy = yes, MRI stage = IA",
           "radiotherapy = yes, MRI stage = IB"),
           xlab="OS Time (days)")
btm<-grid.text("n = 325\n(subset: MRI stage = I)",gp=gpar(fontsize=16))
png("KMsurvCurves_radiot&MRIst.png",width=12,height=8,units='in',res=300)
grid.arrange(gg2$plot,gg2$table,btm,layout_matrix=rbind(c(1,1),c(1,1),c(1,1),c(1,1),c(1,1),c(2,2),c(2,2),c(3,3)))
dev.off()

radioty<-radiot[which(radiot=="yes")]
yess<-which(radiot=="yes")
mainfit.yradiot<-survfit(radiot.surv~radiot+stageRMN,data=dat,subset=yess)
library(survminer)
ggst<-ggsurvplot(mainfit.yradiot,risk.table = TRUE,#pval = TRUE,
           title="Kaplan-Meier estimate with MRI stage unadjusted",#"Before Fitting the Model",
           xlab="OS Time (days)",
           legend.title="MRI stage",
           legend.labs=c("IA", "IB"))
           #font.title=8, font.x=7, font.y=7)
summary(mainfit.yradiot)
#ggst$plot<-ggst$plot+annotate("text", x = 630, y=0.2, label = paste("for log-rank test"),size=5,fontface='bold')
surv_pvalue(
  fit=mainfit.yradiot,
  data = daty,
  method = "survdiff",
  test.for.trend = FALSE,
  combine = TRUE
)
radiot.survy<-radiot.surv[yess]
chemoty<-chemot[yess]
stageRMNy<-stageRMN[yess]
gradehy<-gradeh[yess]
agey<-age[yess]
# md5s<-coxph(radiot.survy~chemoty+stageRMNy+gradehy+agey)
# gst.list
# ggsurvplot_group_by(md5s,data=subset(dat,radiot="yes"),group.by=c("gradeh","stageRMN"))
md55<-coxph(formula = radiot.survy ~ chemoty + stageRMNy + gradehy + agey)
daty<-data.frame(radiot.survy[,1],radiot.survy[,2],chemoty,stageRMNy,gradehy,agey)
fit5y<-surv_fit(md55$formula,data = daty)
p5<-surv_pvalue(
  fit=fit5y,
  data = daty,
  method = "survdiff",
  test.for.trend = FALSE,
  combine = TRUE
)
# surv_pvalue(
#   fit=fit5yy,
#   data = daty[yess,],
#   method = "survdiff",
#   test.for.trend = FALSE,
#   combine = TRUE
# )
p5val<-p5$pval
gga5st<-ggadjustedcurves(md55,data=daty,
                         variable = "stageRMNy",
                         title="Survival curves predicted by the Cox (PH) model with the MRI stage\nadjusted for chemotherapy + tumour grade + age",
                         xlab="OS Time (days)",ylab="Survival probability",
                         method = "average", pval=TRUE,
                         legend="top",legend.title="MRI stage")
gga5st<-gga5st+annotate("text", x = 750, y=0.4, label = paste("p = ",p5val," for log-rank test"),size=5)#,fontface='bold')#+labs=c("A","B")
gga5st<-gga5st+annotate("text", x = 750, y=0.1, label = paste("HR = 1.7859 for MRI stage (IB relative to IA)\n95% CI: 0.85578-3.727"),size=5)
gga5tab<-ggsurvtable(survfit(md55),data=daty, timeby = 1000,xlab="OS Time (days)")
#fit5yy<-survfit(md5,data=dat,subset=yess)
#ggsurvtable(fit5yy,data=dat[yess,], timeby = 1000,xlab="OS Time (days)")#,legend.title = "MRI stage",legend.labs = c("IA","IB"))
summary(md5)
plot(cox.zph(md5))
ggsurv1<-ggsurvplot(survfit(md55),data=dat, ylim=c(.45, 1), xlab="Days since diagnosis", mark.time=F, ylab="Probability of a successful outcome", title="n = 146 (subset: radiotherapy = yes)")#,risk.table = TRUE)
noo<-which(radiot=="no")
radiotn<-radiot[noo]
radiot.survn<-radiot.surv[noo]
chemotn<-chemot[noo]
stageRMNn<-stageRMN[noo]
gradehn<-gradeh[noo]
agen<-age[noo]
md5n<-coxph(radiot.survn ~ chemotn + stageRMNn + gradehn + agen)
ggsurv2<-ggsurvplot(survfit(md5n),data=dat, ylim=c(.45, 1), xlab="Days since diagnosis (OS Time)", mark.time=F, ylab="Probability of a successful outcome", title="n = 179 (subset: radiotherapy = no)")#,risk.table = TRUE)
#family=windowsFont("Arial")),col=c(1, 2, 4)
gglist<-list("radiotherapy = yes (n = 146)"=survfit(md55),"radiotherapy = no (n = 179)"=survfit(md5n))
png("DifferenceBetweenYes&NoRadiot.png",width=16,height=8,units='in',res=300)
gsc<-ggsurvplot_combine(gglist,data=dat,risk.table = TRUE,title = "Did radiotherapy make a difference in patients with MRI stage I according to the prediction of \nthe Cox model with radiotherapy + chemotherapy + MRI stage + tumour grade + age?",
                  ylim=c(.45, 1), xlab="Days since diagnosis (OS Time)", mark.time=F, ylab="Probability of a successful outcome")
grid.arrange(gsc$plot,gsc$table,layout_matrix=rbind(c(1,1),c(1,1),c(1,1),c(1,1),c(2,2)))
dev.off()
lisummary(survfit(md5))
# age.group1<-vector()
# age.group1[which(age<45)]<-'<45'
# age.group1[which(age>=45&age<=70)]<-'45-70'
# age.group1[which(age>70)]<-'>70'
# md5strat<-coxph(radiot.surv ~ strata(radiot,age) + chemot +stageRMN + gradeh)
# ggadjustedcurves(md5strat,data=dat[yess,],
#                  variable = "stageRMN",
#                  title="Cox (PH) model estimate with the MRI stage adjusted for chemotherapy + tumour grade + age",
#                  xlab="OS Time (days)",ylab="Survival probability",
#                  method = "average", pval=TRUE,
#                  legend="top",legend.title="MRI stage",
# )
# surv_diff<- function(df){
#   sdiff <- survdiff(Surv(time, status) ~ sex, data = df)
#   pvalue <- pchisq(sdiff$chisq, length(sdiff$n) - 1, lower.tail = FALSE)
#   return(pvalue)
# }
# gga5<-ggsurvplot(md5,data=dat,
#                          variable = "stageRMN",
#                          title="Cox (PH) model estimate with the MRI stage adjusted for chemotherapy + tumour grade + age",
#                          xlab="OS Time (days)",ylab="Survival probability",
#                          method = "average", pval=TRUE,
#                          legend="top",legend.title="MRI stage",
# )
# gga5$plot<-gga5$plot+geom_point()+facet_grid( .~ stageRMN)

#gga8st<-ggadjustedcurves(md8,data=dat[yess,],variable = "stageRMN",title="Cox (PH) Model with the MRI stage adjusted for chemotherapy ",xlab="Time (days)")
png("compareSurvCurvesGrByMRIstage-Adj&Unadj.png",width=16,height=8,units='in',res=300)
library(grid)
library(gridExtra)
# createBoldTitle <- function(text, gpar, vjust=NULL, hjust=NULL){
#   txt <- textGrob(text, gp=gpar, vjust = vjust, hjust = hjust)
#   #undLine <- segmentsGrob(x0 = txt$x - grobWidth(txt), y0 = txt$y - unit(0.55, "grobheight", txt), x1 = txt$x, y1 = txt$y - unit(0.55, "grobheight", txt))
#   gt <- grobTree(txt, undLine)
# }
library(extrafont)
#loadfonts(device = "win")
#font_import() -to import them (all)
#loadfonts()
fonts()# -to load them
#windowsFonts()

title<-grid.text("The Overall Survival Grouped by the MRI Stage",gp=gpar(fontsize=20,fontface="bold",fontfamily="Calibri"))#,fontfamily="Calibri"))
subtitle<-grid.text("Survival Curves\n",gp=gpar(fontsize=16,fontfamily="Arial"))
margin <- unit(0.5, "line")
#marginb<-unit(-0.5,"line")
#n<-grid.text("n = 146",gp=gpar(fontface="bold"))
#nh<-nchar("n = 146")
#bot<-grid.text(" (subset: radiotherapy = yes, MRI stage = I; IA-IB = unbalanced groups)",gp=gpar(fontsize=14))
#nbot<-nchar(" (subset: radiotherapy = yes, MRI stage = I; IA-IB = unbalanced groups)")
btm<-grid.text("n = 146\n(subset: radiotherapy = yes, MRI stage = I; IA-IB = unbalanced groups)",gp=gpar(fontsize=16))
garr<-grid.arrange(ggst$plot,gga5st,ggst$table,gga5tab$risk.table,layout_matrix = rbind(c(1,2),c(1,2),c(1,2),c(3,4)))
grid.arrange(title,subtitle,garr,
             heights = unit.c(grobHeight(title) + 1.2*margin, 
                              grobHeight(subtitle) + margin, 
                              #grobHeight(bot)+marginb,
                              unit(1,"null")
                              ), bottom=btm)#grid.arrange(n,bot,widths=9:1))#layout_matrix=rbind(c(1,2),c(1,2))))
dev.off()


stageRMN[which(stageRMN=="IA")]
oneA<-which(stageRMN=="IA")
stageRMN[which(stageRMN=="IB")]
oneB<-which(stageRMN=="IB")
mainfit.IAmri<-survfit(radiot.surv~radiot+stageRMN,data=dat,subset=oneA)
title<-grid.text("The Overall Survival Grouped by Radiotherapy\n",gp=gpar(fontsize=20,fontface="bold",fontfamily="Calibri"))
subtitle1<-grid.text("Kaplan-Meier estimates with radiotherapy unadjusted\n",gp=gpar(fontsize=16,fontfamily="Arial"))
ggs1<-ggsurvplot(mainfit.IAmri,risk.table = TRUE,#pval = TRUE,
           xlab="OS Time (days)",
           legend.title="radiotherapy",
           legend.labs=c("yes", "no"))
surv_pvalue(
  fit=mainfit.IAmri,
  data = dat,
  method = "survdiff",
  test.for.trend = FALSE,
  combine = TRUE
)
mainfit.IBmri<-survfit(radiot.surv~radiot+stageRMN,data=dat,subset=oneB)
ggs2<-ggsurvplot(mainfit.IBmri,risk.table = TRUE,#pval = TRUE,
           xlab="OS Time (days)",
           legend.title="radiotherapy",
           legend.labs=c("yes", "no"))
surv_pvalue(
  fit=mainfit.IBmri,
  data = dat,
  method = "survdiff",
  test.for.trend = FALSE,
  combine = TRUE,
  rho=1
)
garr1<-grid.arrange(ggs1$plot,ggs1$table,layout_matrix = rbind(c(1,1),c(1,1),c(1,1),c(2,2)))
btm1<-grid.text("n = 219\n(subset: MRI stage = IA)", gp=gpar(fontsize=16))
garr2<-grid.arrange(ggs2$plot,ggs2$table,layout_matrix = rbind(c(1,1),c(1,1),c(1,1),c(2,2)))
btm2<-grid.text("n = 106\n(subset: MRI stage = IB)", gp=gpar(fontsize=16))
btmf<-grid.arrange(btm1,btm2,ncol=2)#,layout_matrix = rbind(c(1,2)))
grr<-grid.arrange(garr1,garr2,ncol=2)
grrf<-grid.arrange(grr,btmf,layout_matrix=rbind(c(1,1),c(1,1),c(1,1),c(1,1),c(2,2)))
##
subtitle2<-grid.text("Survival curves predicted by the Cox (PH) model with the radiotherapy adjusted for chemotherapy + tumour grade + age",gp=gpar(fontsize=16,fontfamily="Arial"))
radiot.survA=radiot.surv[oneA]
radiotA=radiot[oneA]
chemotA=chemot[oneA]
stageRMNA=stageRMN[oneA]
gradehA=gradeh[oneA]
ageA=age[oneA]
md5A<-coxph(radiot.survA ~ radiotA + chemotA + gradehA + ageA)
datA<-data.frame(radiot.survA[,1],radiot.survA[,2],radiotA,chemotA,stageRMNA,gradehA,ageA)
gga5rad<-ggadjustedcurves(md5A,data=datA,
                         variable = "radiotA",
                         xlab="OS Time (days)",ylab="Survival probability",
                         method = "average", pval=TRUE,
                         legend="top",legend.title="radiotherapy")
fit5A<-surv_fit(md5A$formula,data = datA)
p5<-surv_pvalue(
  fit=fit5A,
  data = datA,
  method = "survdiff",
  test.for.trend = FALSE,
  combine = TRUE
)
p5val<-p5$pval
gga5rad<-gga5rad+annotate("text", x = 750, y=0.4, label = paste("p = ",p5val," for log-rank test"),size=5)#,fontface='bold')#+labs=c("A","B")
gga5rad<-gga5rad+annotate("text", x = 750, y=0.1, label = paste("HR = 0.1339 for radiotherapy (yes relative to no)\n95% CI: 0.05365-0.3341"),size=5)
gga5tab<-ggsurvtable(survfit(md5A),data=datA, timeby = 1000,xlab="OS Time (days)")
##
radiot.survB=radiot.surv[oneB]
radiotB=radiot[oneB]
chemotB=chemot[oneB]
stageRMNB=stageRMN[oneB]
gradehB=gradeh[oneB]
ageB=age[oneB]
md5B<-coxph(radiot.survB ~ radiotB + chemotB + gradehB + ageB)
datB<-data.frame(radiot.survB[,1],radiot.survB[,2],radiotB,chemotB,stageRMNB,gradehB,ageB)
gga5radB<-ggadjustedcurves(md5B,data=datB,
                          variable = "radiotB",
                          xlab="OS Time (days)",ylab="Survival probability",
                          method = "average", pval=TRUE,
                          legend="top",legend.title="radiotherapy")
fit5B<-surv_fit(md5B$formula,data = datB)
p5<-surv_pvalue(
  fit=fit5B,
  data = datB,
  method = "survdiff",
  test.for.trend = FALSE,
  combine = TRUE
)
p5val<-p5$pval
gga5radB<-gga5radB+annotate("text", x = 750, y=0.4, label = paste("p = ",p5val," for log-rank test"),size=5)#,fontface='bold')#+labs=c("A","B")
gga5radB<-gga5radB+annotate("text", x = 750, y=0.1, label = paste("HR = 0.3425 for radiotherapy (yes relative to no)\n95% CI: 0.1118-1.049"),size=5)
gga5tabB<-ggsurvtable(survfit(md5B),data=datB, timeby = 1000,xlab="OS Time (days)")

gr5<-grid.arrange(top=subtitle2,gga5rad,gga5radB,gga5tab$risk.table,gga5tabB$risk.table,layout_matrix=rbind(c(1,2),c(1,2),c(1,2),c(1,2),c(3,4)))#(c(1,1),c(2,3),c(2,3),c(2,3),c(2,3),c(4,5)))#(c(1,1),c(2,2),c(2,2),c(2,2),c(2,2),c(2,2),

subtitle3<-grid.text("\nSurvival curves predicted by the Cox (PH) model with the radiotherapy adjusted for chemotherapy + CRS + age",gp=gpar(fontsize=16,fontfamily="Arial"))
CRSA<-CRS[oneA]
md8A<-coxph(radiot.survA ~ radiotA + chemotA + CRSA + ageA)
dat8A<-data.frame(radiot.survA[,1],radiot.survA[,2],radiotA,chemotA,stageRMNA,CRSA,ageA)
gga8rad<-ggadjustedcurves(md8A,data=dat8A,
                          variable = "radiotA",
                          xlab="OS Time (days)",ylab="Survival probability",
                          method = "average", pval=TRUE,
                          legend="top",legend.title="radiotherapy")
fit8A<-surv_fit(md8A$formula,data = dat8A)
p8<-surv_pvalue(
  fit=fit8A,
  data = dat8A,
  method = "survdiff",
  test.for.trend = FALSE,
  combine = TRUE
)
p8val<-p8$pval
gga8rad<-gga8rad+annotate("text", x = 750, y=0.4, label = paste("p = ",p8val," for log-rank test"),size=5)#,fontface='bold')#+labs=c("A","B")
gga8rad<-gga8rad+annotate("text", x = 750, y=0.1, label = paste("HR = 0.2499 for radiotherapy (yes relative to no)\n95% CI: 0.08877-0.7036"),size=5)
gga8tab<-ggsurvtable(survfit(md8A),data=dat8A, timeby = 1000,xlab="OS Time (days)")
##
CRSB<-CRS[oneB]
md8B<-coxph(radiot.survB ~ radiotB + chemotB + CRSB + ageB)
dat8B<-data.frame(radiot.survB[,1],radiot.survB[,2],radiotB,chemotB,stageRMNB,CRSB,ageB)
gga8radB<-ggadjustedcurves(md8B,data=dat8B,
                          variable = "radiotB",
                          xlab="OS Time (days)",ylab="Survival probability",
                          method = "average", pval=TRUE,
                          legend="top",legend.title="radiotherapy")
fit8B<-surv_fit(md8B$formula,data = dat8B)
p8<-surv_pvalue(
  fit=fit8B,
  data = dat8B,
  method = "survdiff",
  test.for.trend = FALSE,
  combine = TRUE
)
p8val<-p8$pval
gga8radB<-gga8radB+annotate("text", x = 750, y=0.4, label = paste("p = ",p8val," for log-rank test"),size=5)#,fontface='bold')#+labs=c("A","B")
gga8radB<-gga8radB+annotate("text", x = 750, y=0.1, label = paste("HR = 0.24371 for radiotherapy (yes relative to no)\n95% CI: 0.076605-0.7753"),size=5)
gga8tabB<-ggsurvtable(survfit(md8B),data=dat8B, timeby = 1000,xlab="OS Time (days)")        

gr8<-grid.arrange(top=subtitle3,gga8rad,gga8radB,gga8tab$risk.table,gga8tabB$risk.table,layout_matrix=rbind(c(1,2),c(1,2),c(1,2),c(1,2),c(3,4)))                             

gr58<-grid.arrange(gr5,gr8,nrow=2)
#grf<-grid.arrange(grrf,gr58,layout_matrix=rbind(rep(c(1,1),8),rep(c(2,2),16)))
grf<-grid.arrange(gr58,btmf,layout_matrix=rbind(c(1,1),c(1,1),c(1,1),c(1,1),c(1,1),c(1,1),c(1,1),c(1,1),c(1,1),c(2,2)))
                                                                   #c(3,4),c(3,4),c(3,4),c(3,4),c(5,6)))
png("stageIAvsIB-w&woutRadiotModels.png",width=18,height=15,units='in',res=300)
grid.arrange(title,grf,
             heights = unit.c(grobHeight(title) + 1.2*margin, 
                              #grobHeight(subtitle) + margin, 
                              unit(1,"null")))
dev.off()

fit5_8<-list("model 5 fit"=surv_fit(md5$formula,dat),"model 8 fit"=surv_fit(md8$formula,dat))
#it also works with survfit function (i.e. without underscore)
spc<-surv_pvalue(
  fit=fit5_8,
  data = dat,
  method = "survdiff",
  test.for.trend = FALSE,
  combine = TRUE
)
pvalc<-spc$pval
fit5_8_<-list("model 5 fit"=survfit(md5),"model 8 fit"=survfit(md8))
ggc<-ggsurvplot_combine(fit5_8_,dat,title="Comparing the Predicted Survival Functions for the Two Cox (PH) Models",xlab="OS Time (days)",pval=TRUE,pval.method=TRUE)#plotted together

fit5<-surv_fit(md5$formula,data = dat)
sp.md5<-surv_pvalue(
  fit=fit5,
  method = "survdiff",
  combine = FALSE
)
pval.md5<-sp.md5$pval
survdiff(md5$formula,data = dat)
gga5<-ggadjustedcurves(md5,data=dat,
                 title="The Predicted Survival Function for Model 5 with\nradiotherapy + chemotherapy + MRI stage + tumour grade + age",
                 xlab="OS Time (days)",
                 ylab="Survival probability")#now plotted separately 
gga5<-gga5+annotate("text", x = 500, y=0.05, label = paste("p = ",pval.md5),size=5,fontface='bold')
gga5tab<-ggsurvtable(survfit(md5),data=dat, timeby = 1000,xlab="OS Time (days)")

fit8<-surv_fit(md8$formula,data = dat)
sp.md8<-surv_pvalue(
  fit=fit8,
  method = "survdiff",
  combine = FALSE
)
pval.md8<-sp.md8$pval
survdiff(md8$formula,data = dat)
gga8<-ggadjustedcurves(md8,data=dat,
                 title="The Predicted Survival Function for Model 8 with radiotherapy + chemotherapy + CRS + age",
                 xlab="OS Time (days)",ylab="Survival probability")
gga8<-gga8+annotate("text", x = 500, y=0.05, label = paste("p = ",pvalc[2]),size=5,fontface='bold')
gga8tab<-ggsurvtable(survfit(md8),data=dat, timeby = 1000,xlab="OS Time (days)")
btm<-grid.text("n = 325\n(subset: MRI stage = I)",gp=gpar(fontsize=16))
png("compareSurvF-betweenModels-squaredGrid.png",width=18,height=12,units='in',res=300)
grid.arrange(ggc$plot,gga5,gga8,gga5tab$risk.table,gga8tab$risk.table,layout_matrix = rbind(c(1,1),c(1,1),c(1,1),c(2,3),c(2,3),c(2,3),c(4,5)),
             bottom=btm)
dev.off()

ggadjustedcurves(md5,data=dat,variable = "radiot")
ggadjustedcurves(md8,data=dat,variable = "radiot")
#our model tells us though that survival increases with radiotherapy
basehaz(md8)

#risk.table = FALSE,pval = TRUE
# radiot1<-radiot
# radiot1<-factor(radiot1)#or could do radiot1<-factor(radiot) directly
# chemot1<-chemot
# chemot1<-factor(chemot1)
# stageRMN1<-stageRMN
# stageRMN1<-factor(stageRMN1)
# gradeh1<-gradeh
# gradeh1<-factor(gradeh1)
# CRS1<-CRS
# CRS1<-factor(CRS1)
# dat1<-data.frame(stime,sevent,radiot.surv,radiot1,chemot1,stageRMN1,gradeh1,age,CRS1)
# md5f<-coxph(radiot.surv ~ radiot1 + chemot1 + gradeh1 + age)
# md5y<-survfit(md5f,subset = yess)
# colnames(dat1)[6]=factor(colnames(dat1)[6])
# vv=vector()
# vv<-c("stageRMN1","chemot1")
# ggsurvplot_group_by(md5y,data=dat1,group.by ="stageRMN1")
