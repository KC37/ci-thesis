KM.crs<-function(){
#which(grepl("\\?",data$Death.Date)==FALSE))
  #op.date<-data$OpDate[inclInStudy]
  #death.date<-data$Death.Date[inclInStudy]
  # capture the characters after the first space
  #sub(".*? ", "", data$Death.Date[inclInStudy])
  
  #gsub("(.+?)(\\_.*)", "\\1", data$Death.Date[inclInStudy])
  #capture the characters before first space 
  #(star position against the point means before or after)
  ddsub<-sub(" .*", "", data$Death.Date[inclInStudy])
  #sub(" .*", "","no ") returns "no"
  odate<-sub(" .*", "", data$OpDate[inclInStudy])
  # #sub(" .*", "", "2018-01-25 00:00:00 hey")
  # library("sjmisc")
  ddate<-ifelse(ddsub=="no" | ddsub=="No", "2020-08-03", ddsub)
  # table(death.date)
  # ddvalid<-grepl("-",death.date)
  # ddate<-death.date[ddvalid]
  # odate<-op.date[ddvalid]
  
  os.time<-difftime(ddate, odate, tz="GMT", units = "days")
  os.event<-as.numeric(ddate!="2020-08-03")#is.deceased
  
  #see if there are other dates in the same year as the end of the study
  #ddate[which(grepl("2020",ddate)&ddate!="2020-08-03")]
  risk.score<-data$Risk.score[setdiff(inclInStudy,which(data$Risk.score=="x"))]#which(data$Reasons.for.Exclusion=="include"&(grepl("-",data$Death.Date)|data$Death.Date %in% c("no","No","no ")))]
  data$Death.Date[intersect(inclInStudy,which(data$Risk.score=="x"))]
  data$Reasons.for.Exclusion[intersect(inclInStudy,which(data$Risk.score=="x"))]
  #159
  risks<-data$Risk.score[inclInStudy]
  i<-which(risks=="x")
  os.time<-os.time[-i]
  os.event<-os.event[-i]
  
  library(survival)
  ec.os<-Surv(os.time,os.event)
  
  # ddd<-data$Death.Date[which(data$`Reasons for Exclusion`=="include"&(grepl("-",data$Death.Date)|data$Death.Date %in% c("no","No")))]
  # ddd.fin<-ifelse(ddd=="no" | ddd=="No","2020-08-03",ddd)
  # ddd.fins<-sub(" .*", "", ddd.fin)
  # setdiff(ddate,ddd.fins)
  # for(i in 1:length(ddate))
  #   if(ddate[i]!=ddd.fins[i])
  #     print(i)
  # data$Death.Date[ddate!=ddd.fins]
  length(risk.score)
  length(inclInStudy)
  length(which(data$Reasons.for.Exclusion=="include"))
  risk.score<-ifelse(risk.score=="Intermediate","intermediate",risk.score)
  # should do that to the whole table actually
  #data$Risk.score<-ifelse(data$Risk.score=="Intermediate","intermediate",data$Risk.score)
  #length(data$Risk.score)
  length(which(risk.score=="low"))
  length(which(risk.score=="intermediate"))
  length(which(risk.score=="high"))
  length(which(risk.score=="advanced"))
  rsc<-factor(risk.score,levels=c("low","intermediate","high","advanced"))
  png("ggplot_EC_OS_byRiskScore304.png",width=16,height=8,units='in',res=300)
  library(survminer)
  #library(RColorBrewer)
  ggsurvplot(
    survfit(ec.os ~ rsc),
    title = "The Overall Survival Based on Clinical Risk Score",
    #subtitle = "Stratification by Risk Score",
    data = data,
    legend.title = "CR Score",
    risk.table = TRUE,
    pval = TRUE,
    conf.int = FALSE, # show confidence intervals for 
                      # point estimates of survival curves.
    xlim = c(0,5000),
    xlab = "OS Time (days)",
    break.time.by = 500,
    #ggtheme = theme_light(),
    #brewer.pal(4,"Set2")[3]
    palette = rev(c("red","orange","deepskyblue","forestgreen")),
    risk.table.title="Strata (by CR Score) - size every 500 days",
    #risk.table.ticks.col = TRUE,
    #risk.table.y.text.col = T,
    #risk.table.y.text = FALSE,
    legend.labs = levels(rsc)#c("advanced","high","intermediate","low")
    #legend.labs = c("risk.score=advanced", "risk.score=high", "risk.score=intermediate", "risk.score=low")
  )
  dev.off()
  
  min(os.time)
  max(os.time)
  
  surv.risk.table<-data.frame(inclInStudy[-i],risk.score,os.time,os.event)
  surv.risk.table[order(surv.risk.table$os.time,decreasing=TRUE),]
  #the next advanced risk score has shorter survival that at least one patient of each/every other type
  c(data$PatientID[2],data$OpDate[2],data$Death.Date[2],data$Risk.score[2])
  #1 outlier: joined the study in 2007 and was censored
  
  
  # low.rs<-which(risk.score=="low")
  # adv.rs<-which(risk.score=="advanced")
  # max(os.time[low.rs])
  # max(os.time[adv.rs])
  # high.risk <- ifelse(risk.score=="high" | risk.score=="advanced",1,0)
  # #as.numeric(risk.score>median(risk.score))
  # png("ggplot_EC_OS_by2Risk.png",width=16,height=8,units='in',res=300)
  # library(survminer)
  # ggsurvplot(
  #   survfit(ec.os ~ high.risk),
  #   data = data,
  #   risk.table = TRUE,
  #   pval = TRUE,
  #   conf.int = TRUE,
  #   xlim = c(0,5000),
  #   xlab = "OS Time (days)",
  #   break.time.by = 500,
  #   ggtheme = theme_light(),
  #   risk.table.y.text.col = T,
  #   risk.table.y.text = FALSE
  #   #legend.labs = c("risk.score=advanced", "risk.score=lo", "risk.score=high", "risk.score=advanced")
  # )
  # dev.off()
  }
  
  
  
