age.period<-function(){
  sort(data$AgeAtDiagnosis[allincl])#27-93
}


determine.time.period<-function(){
#long-term survival
dd<-ddate1[-213]
sort(as.Date(dd, format="%Y-%m-%d"))
sort(as.Date(opd, format="%Y-%m-%d"))
#2013-2020 for death
#2007-2018 for diagnostic (surgery/operation date)
#plot created before gathering the final study population
png(file = "nysurv-barchart271.png",width=10,height=8,units='in',res=300)
one<-n.year.surv(1)
three<-n.year.surv(3)
five<-n.year.surv(5)
ten<-n.year.surv(10)
H <- c(one,three,five,ten)
M <- c("One","Three","Five","Ten")
barplot(H,names.arg=M,
      main="Endometrial Cancer One-, Three-, Five- and Ten-Year Net Survival, Adults (Aged 27-93), 2007-2018",
      xlab="Years after Diagnosis",ylab="Net Survival (%)",
      #xlim=c(0,100),
      col="magenta"
     )
dev.off()

#after establishing the final study population
#plot generated with pre-built function in R
png(file = "nysurv-barchart458.png",width=16,height=8,units='in',res=300)
one<-n.year.surv(1,allstudypop)
three<-n.year.surv(3,allstudypop)
five<-n.year.surv(5,allstudypop)
#censored;ten<-n.year.surv(10,allstudypop)
H <- c(one,three,five)
M <- c("One","Three","Five")
barplot(H,names.arg=M,
        main="Endometrial Cancer One-, Three-, and Five-Year Survival, Adults (Aged 27-97 and Diagnosed during 2008-2018)",
        xlab="Years after Diagnosis",ylab="Survival Rate (%)",
        ylim=c(0,100),
        col="magenta"
)
dev.off()

#same plot generated with ggplot2
library(ggplot2)
nysurv.data<-data.frame(H,M)
M1<-M
M<-factor(M)
nysurv.data$M <- factor(nysurv.data$M, levels = c("One","Three","Five"))
png(file = "nysurv-ggplot-barchart458.png",width=16,height=8,units='in',res=300)
ggplot(nysurv.data,aes(x = M, y = H)) +
  geom_col(aes(fill = "magenta")) +
  geom_label(aes(label = round(H, 2))) +
  xlab("Years after Diagnosis") +
  ylab("Survival Rate (%)")+ylim(c(0,100))+
  ggtitle("Endometrial Cancer One-, Three-, and Five-Year Survival, Adults (458 Patients Aged 27-97 and Diagnosed during 2008-2018)")+
  theme(legend.position = "none")
dev.off()

surv.obj<-generate.Surv.obj(allstudypop)
#estimating 1, 3, 5, and 10-year survival using Survival object
#incorrect estimate of the 1-year probability of survival when ignoring the fact that some patients were censored before 1 year
oney<-summary(survfit(surv.obj ~ 1, data = corrected.data), times = 365.25)
threey<-summary(survfit(surv.obj ~ 1, data = corrected.data), times = 3*365.25)
fivey<-summary(survfit(surv.obj ~ 1, data = corrected.data), times = 5*365.25)
all.data<-corrected.data[allstudypop,]
summary(survfit(surv.obj ~ 1, data = data.frame(surv.obj)), times = 10*365.25)
oney
#365    434      24    0.948  0.0104        0.927        0.968
summary(survfit(surv.obj ~ 1, data = all.data), times = 365.25)

#plot generated with ggplot2 with the n-year survival calculated by R, which is similar to the results of my computation (plotted above)
ysurv<-c(oney$surv*100,threey$surv*100,fivey$surv*100)
nams<-M[-4]
nysurv.dat<-data.frame(ysurv,nams)
nysurv.dat$nams <- factor(nysurv.dat$nams, levels = c("One","Three","Five"))
png(file = "nysurv-barchart458.png",width=16,height=8,units='in',res=300)
ggplot(nysurv.dat,aes(x = nams, y = ysurv)) +
  geom_col(aes(fill = "magenta")) +
  geom_label(aes(label = round(ysurv, 2))) +
  xlab("Years after Diagnosis") +
  ylab("Survival Rate (%)")+ylim(c(0,100))+
  ggtitle("Endometrial Cancer One-, Three-, and Five-Year Survival (458 Patients Aged 27-97 and Diagnosed during 2008-2018)")+
  theme(legend.position = "none")
dev.off()
}
