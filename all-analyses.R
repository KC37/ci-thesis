# the four analyses developed for the Cancer Informatics Thesis

# comparison between stages IA (superficial) and IB (deep)
#incorrectly staging based on the depth of myometrial invasion

#loading the necessary libraries
library(sjmisc) # for str_contains
library(qmrparser) # for isDigit
library(useful) # for upper.case
library(stringr) # for str_replace, str_extract

# initializing the variables we will be using (cols+df)
init.vars<-function(){

  pid<-data$PatientID
  pid

  histo.stage<-data$STAGE
  table(histo.stage)
  
  # already merged MRI staging data (int+ext)
  #found under FINAL STAGE, MRI stage
  MRI.stage<-data$"X.2"
  table(MRI.stage)

  # merging internal and external MRI depth
  MRI.minv <- ifelse(data$"Sup.Deep"!="x" & data$"Sup.Deep"!="", data$"Sup.Deep", data$"Sup.Deep.1")
  table(MRI.minv)
  
  histo.minv<-data$"Sup.Deep.2"
  table(histo.minv)
  # this was for stage I: IA vs IB
  #library(h2o)
  stIAvalid<-which(data$X.2=="IA"&str_contains(data$Sup.Deep,"sup",ignore.case = TRUE)|
                     data$X.2=="IA"&str_contains(data$Sup.Deep.1,"sup",ignore.case = TRUE)&!str_contains(data$Sup.Deep.1,"deep",ignore.case = TRUE) )
  cat(length(stIAvalid),"out of",length(which(data$X.2=="IA")),"cases are valid, i.e. the MRI findings are in accordance with the MRI staging (IA)",sep=" ")
  
  stIBvalid<-which(data$X.2=="IB"&str_contains(data$Sup.Deep,"deep",ignore.case = TRUE)|
                     data$X.2=="IB"&str_contains(data$Sup.Deep.1,"deep",ignore.case = TRUE)&!str_contains(data$Sup.Deep.1,"sup",ignore.case = TRUE) )
  cat(length(stIBvalid),"out of",length(which(data$X.2=="IB")),"cases are valid, i.e. the MRI findings are in accordance with the MRI staging (IB)",sep=" ")
  
  MRI.strinv <- data$Y.N.4
  table(MRI.strinv)
  histo.strinv <- data$Y.N.5
  table(histo.strinv)
  # this was for stage II
  stIIvalid<-which(data$X.2=="II"&tolower(data$Y.N.4)=="y")
  cat(length(stIIvalid),"out of",length(which(data$X.2=="II")),"cases are valid, i.e. the MRI findings are in accordance with the the MRI staging (II)",sep=" ")
  
  MRI.vagpar.met <- data$Y.N.6
  table(MRI.vagpar.met)
  histo.vagpar.met <- data$Y.N.7
  table(histo.vagpar.met)
  # this was for stage IIIB
  stIIIBvalid<-which(data$X.2=="IIIB"&tolower(data$Y.N.6)=="y")
  discrepIIIB<-which(data$X.2=="IIIB"&tolower(data$Y.N.6)=="n")
  cat(length(stIIIBvalid),"out of",length(stIIIBvalid)+length(discrepIIIB),"cases are valid, i.e. the MRI findings are in accordance with the MRI staging (IIIB)",sep=" ")
  #do not show discrepancies between MRI observations/findings and staging (IIIB)
  data$Y.N.6[which(data$X.2=="IIIB")]
  missing.dataIIIB<-which(data$X.2=="IIIB"&tolower(data$Y.N.6)!="y"&tolower(data$Y.N.6)!="n")
  length(which(data$X.2=="IIIB"))
  #no misssing data for stage IIIB
  
  
  MRI.adinv <- data$Y.N.8
  table(MRI.adinv)
  histo.adinv <- data$Y.N.9
  table(histo.adinv)
  MRI.serinv <- data$Y.N.10
  table(MRI.serinv)
  histo.serinv <- data$Y.N.11
  table(histo.serinv)
  # these were for stage IIIA
  stIIIAvalid<-which(data$X.2=="IIIA"&tolower(data$Y.N.8)=="y"|data$X.2=="IIIA"&tolower(data$Y.N.10)=="y")
  discrepIIIA<-which(data$X.2=="IIIA"&tolower(data$Y.N.8)=="n"&tolower(data$Y.N.10)=="n")
  cat(length(stIIIAvalid),"out of",length(stIIIAvalid)+length(discrepIIIA),"cases are valid, i.e. the MRI findings are in accordance with the MRI staging (IIIA)",sep=" ")
  data$Y.N.8[which(data$X.2=="IIIA")]
  data$Y.N.10[which(data$X.2=="IIIA")]
  missing.dataIIIA<-which(data$X.2=="IIIA"&tolower(data$Y.N.8)!="y"&tolower(data$Y.N.10)!="y"&tolower(data$Y.N.8)!="n"&tolower(data$Y.N.10)!="n")
  length(which(data$X.2=="IIIA"))
  #no misssing data for stage IIIA
  
  
  MRI.peritinv <- data$Y.N.12
  table(MRI.peritinv)
  histo.peritinv <- data$Y.N.13
  table(histo.peritinv)
  # this was for stage IV
  stIVvalid<-which(data$X.2=="IV"&tolower(data$Y.N.12)=="y")
  discrepIV<-which(data$X.2=="IV"&tolower(data$Y.N.12)=="n")
  #the rest is missing data
  cat(length(stIVvalid),"out of",length(stIVvalid)+length(discrepIV),"cases are valid, i.e. the MRI findings are in accordance with the MRI staging (IV)",sep=" ")
  data$Y.N.12[which(data$X.2=="IV")] #len=10
  which(data$X.2=="IV") #len 10
  data$PatientID[which(data$X.2=="IV"&tolower(data$Y.N.12)=="n")]
  data$MRI.date[which(data$X.2=="IV"&tolower(data$Y.N.12)=="n")]
  data$AgeAtDiagnosis[which(data$X.2=="IV"&tolower(data$Y.N.12)=="n")]
  missing.dataIV<-which(data$X.2=="IV"&tolower(data$Y.N.12)!="y"&tolower(data$Y.N.12)!="n")
  data$Y.N.12[missing.dataIV]
  data$X.2[missing.dataIV]
  
  
  MRI.RpelvLN <- data$Y.N.14
  table(MRI.RpelvLN)
  histo.RpelvLN <- data$Y.N.15
  table(histo.RpelvLN)
  MRI.LpelvLN <- data$Y.N.16
  table(MRI.LpelvLN)
  histo.LpelvLN <- data$Y.N.17
  table(histo.LpelvLN)
  # these were for stage IIIC1
  stIIIC1valid<-which(data$X.2=="IIIC1"&tolower(data$Y.N.14)=="y"|data$X.2=="IIIC1"&tolower(data$Y.N.16)=="y")
  discrepIIIC1<-which(data$X.2=="IIIC1"&tolower(data$Y.N.14)=="n"&tolower(data$Y.N.16)=="n")
  cat(length(stIIIC1valid),"out of",length(which(data$X.2=="IIIC1")),"cases are valid, i.e. do not show discrepancies between MRI findings and staging (IIIC1)",sep=" ")
  
  
  MRI.paraortLN <- data$Y.N.18
  table(MRI.paraortLN)
  histo.paraortLN <- data$Y.N.19
  table(histo.paraortLN)
  # these were for stage IIIC2
  stIIIC2valid<-which(data$X.2=="IIIC2"&tolower(data$Y.N.18)=="y")
  discrepIIIC2<-which(data$X.2=="IIIC2"&tolower(data$Y.N.18)=="n")
  cat(length(stIIIC2valid),"out of",length(which(data$X.2=="IIIC2")),"cases are valid, i.e. do not show discrepancies between MRI findings and staging (IIIC2)",sep=" ")
  
  
  MRI.RingnLN <- data$Y.N.20
  table(MRI.RingnLN)
  histo.RingnLN <- data$Y.N.21
  table(histo.RingnLN)
  MRI.LingnLN <- data$Y.N.22
  table(MRI.LingnLN)
  histo.LingnLN <- data$Y.N.23
  table(histo.LingnLN)
  # these were for stage IVB
  stIVBvalid<-which(data$X.2=="IVB"&tolower(data$Y.N.20)=="y"|data$X.2=="IVB"&tolower(data$Y.N.22)=="y"|data$X.2=="IVB"&tolower(data$Y.N.12)=="y")
  discrepIVB<-which(data$X.2=="IVB"&tolower(data$Y.N.20)=="n"&tolower(data$Y.N.22)=="n"&tolower(data$Y.N.12)=="n")
  cat(length(stIVBvalid),"out of",length(stIVBvalid)+length(discrepIVB),"cases are valid, i.e. do not show discrepancies between MRI findings and staging (IVB)",sep=" ")
  data$Y.N.20[which(data$X.2=="IVB")]
  data$Y.N.22[which(data$X.2=="IVB")]
  cat(length(discrepIVB),"out of",length(stIVBvalid)+length(discrepIVB),"cases show discrepancies between MRI findings and staging (IVB)",sep=" ")
  data$PatientID[which(data$X.2=="IVB"&tolower(data$Y.N.20)=="n"&tolower(data$Y.N.22)=="n")]
  length(which(data$X.2=="IVB"))
  missing.dataIVB<-which(data$X.2=="IVB"&tolower(data$Y.N.20)!="y"&tolower(data$Y.N.20)!="n"&tolower(data$Y.N.22)!="y"&tolower(data$Y.N.22)!="n"&!(tolower(data$Y.N.12) %in% c("y","n")))
  data$Y.N.20[missing.dataIVB]
  data$Y.N.22[missing.dataIVB]
  data$X.2[missing.dataIVB]
  
  
  indices<-c(which(data$X.2=="IV"&tolower(data$Y.N.12)=="n"),which(data$X.2=="IVB"&tolower(data$Y.N.20)=="n"&tolower(data$Y.N.22)=="n"))
  indices
  PatientID<-data$PatientID[indices]
  MRI.date<-data$MRI.date[indices]
  AgeAtDiagnosis<-data$AgeAtDiagnosis[indices]
  MRI.Stage<-data$X.2[indices]
  MRI.myom.inv<-MRI.minv[indices]
  dfrr<-data.frame(PatientID,MRI.date,AgeAtDiagnosis,MRI.Stage)
  dfrr$"Myometrial invasion"<-MRI.myom.inv
  dfrr
  start<-grep("Y.N.4", colnames(data))
  end<-grep("Y.N.22", colnames(data))
  cn<-c("Other involvement: Cervical stroma","Vagina/Parametria","Adnexa","Serosa","Peritoneum/Omentum","R-pelvic LN","L-pelvic LN","Para-aortic LN","R-inguinal LN","L-inguinal LN")
  length(cn)
  dfrr[cn]<-data[indices,seq(start,end,by=5)]
  dfrr
  first<-which(data$X.2=="IV"&tolower(data$Y.N.12)=="n")
  last<-which(data$X.2=="IVB"&tolower(data$Y.N.20)=="n"&tolower(data$Y.N.22)=="n")
  stIV<-"peritoneum/omentum"#cn[5]
  stIVB<-"inguinal lymph nodes"
  first.err<-paste("stage",data$X.2[first],"with no",stIV,"involvement",sep=" ")
  last.err<-paste("stage",data$X.2[last],"with no",stIVB,"involvement",sep=" ")
  dfrr$Error<-c(first.err,last.err)
  dfrr
  data[77,seq(start,end,by=5)]
  
  library(writexl)
  write_xlsx(dfrr,"C:\\Users\\Casi\\Documents\\EndoProj\\MRIdiscrepancies.xlsx")
  
  
  grep("Y.N.4", colnames(data))
  data[,74]
  which(colnames(data) == "Y.N.4")
  
  #library(writexl)
  #write_xlsx(dfm,"C:\\Users\\Casi\\Documents\\EndoProj\\MissingMRIfindingsButStaged.xlsx")
  # insert rows at the end of Excel table
  library(xlsx)
  new.workbook <- loadWorkbook("MRIdiscrepancies.xlsx")
  new.sheets <- getSheets(new.workbook)
  new.sheets
  missing.data<-c(missing.dataIV,missing.dataIVB)
  missing.data
  f<-which(data[missing.data,]$X.2=="IV")#missing.data[1]
  l<-which(data[missing.data,]$X.2=="IVB")#missing.data[2]
  ferr<-paste("staged",data$X.2[f],"even if data about",stIV,"involvement is missing",sep=" ")
  lerr<-paste("staged",data$X.2[l],"even if data about",stIVB,"involvement is missing",sep=" ")
  dfm<-data.frame(data$PatientID[missing.data],data$MRI.date[missing.data],data$AgeAtDiagnosis[missing.data],data$X.2[missing.data],MRI.minv[missing.data],data[missing.data,seq(start,end,by=5)],c(ferr,lerr))
  class(dfm)
  ncol(dfrr)
  ncol(dfm)
  dfrr
  colnames(dfm)<-colnames(dfrr)
  dfall<-rbind(dfrr,dfm)
  dfall
  
  addDataFrame(dfm, new.sheets$Sheet1, startColumn = 1, startRow = 10, row.names=FALSE,col.names = FALSE)
  saveWorkbook(new.workbook, 'MRI_discrepancies.xlsx')
  
  accuracy.data<-data.frame(pid,histo.stage,MRI.stage,histo.minv,MRI.minv)
  head(accuracy.data)
  
  allind<-c(indices,missing.data)
  allind
  rid.df<-data.frame(Radiomics.CaseID=data$radiomicsCaseID[allind])
  rid.df
  insert.col.ExcelTable(rid.df,ncol(dfall)+1,"MRI_discrepancies.xlsx","MRI_all_discrepancies.xlsx")
}

insert.col.ExcelTable<-function(df,scol,inFile,outFile){
  library(xlsx)
  # insert column at the end of Excel table
  nworkbook <- loadWorkbook(inFile)
  nsheets <- getSheets(nworkbook)
  nsheets
  addDataFrame(df, nsheets$Sheet1, startColumn = scol, startRow = 2, row.names=FALSE)
  saveWorkbook(nworkbook,outFile)
}

# cleaning the staging data
clean.stage<-function(){
  i=0
  vec=vector()
  # do not use histo.stage alone bc it remains the same, while accuracy.data does not
  for(stage in accuracy.data$histo.stage){
    i=i+1
    sstage<-toString(stage)
    # introduce the character as it is - in str_contains, 
    #but escape it in str_replace(_all)
    if(str_contains(sstage," ")){
      sstage<-str_replace(sstage,"\\s.+","")
    }
    if(str_contains(sstage,"(")){
      sstage<-str_replace(sstage,"\\s*\\([^\\)]+\\)","")
    }
    if(str_contains(sstage,"?")){ # do not have to escape
      sstage<-str_replace_all(sstage,"\\?","") # escape it
    }
    if(str_contains(sstage,",")){
      sstage<-str_replace(sstage,",","")
    }
    if(sstage=="" || str_contains(sstage,"n",ignore.case=TRUE) || sstage=='x'){
      vec<-c(vec,i)
      print(sstage)
      next # instead of continue
    }
    first.chr=substring(sstage,1,1)
    last.chr=substring(sstage,nchar(sstage),nchar(sstage))
    if(isDigit(last.chr))
      if(nchar(sstage)==1)
        ante.chr<-last.chr
    else{
      ante.chr<-first.chr
      last.chr<-substring(sstage,nchar(sstage)-1,nchar(sstage)-1)
    }
    else
      ante.chr=substring(sstage,nchar(sstage)-1,nchar(sstage)-1)
    if(isDigit(ante.chr)){
      nb<-as.numeric(ante.chr)
      if(nb<=3){
        repl<-paste(replicate(nb, "I"), collapse = "")
      } else if(nb==4){
        repl<-"IV"
      } else{
        repl<-""
      }
      class(ante.chr)
      sstage<-str_replace(sstage,ante.chr,repl)
    }
    if(!isDigit(last.chr)&&!upper.case(last.chr)){
      bigl<-toupper(last.chr)
      sstage<-str_replace(sstage,last.chr,bigl)
    }
    if(sstage!=toString(stage)){
      accuracy.data$histo.stage[i]<-sstage
    }
  }
  print(i)
  vec
  # removing uncleaned histological stages
  accuracy.data<-accuracy.data[-vec,]
  table(accuracy.data$histo.stage)
  table(accuracy.data$MRI.stage) # 1 empty cell missed
  nrow(accuracy.data)
  # removing rows with empty MRI stages
  accuracy.data<-subset(accuracy.data,MRI.stage!="")
  table(accuracy.data$MRI.stage)
}

# cleaning the stage-cleaned data on depth of myometrial invasion
clean.minv<-function(){
  library(Xmisc)
  
  table(accuracy.data$histo.minv)
  accuracy.data$histo.minv[startswith(accuracy.data$histo.minv,"sup",ignore.case=TRUE)]<-"superficial"
  accuracy.data$histo.minv[startswith(accuracy.data$histo.minv,"deep",ignore.case=TRUE)]<-"deep"
  table(accuracy.data$histo.minv)
  
  table(accuracy.data$MRI.minv)
  #be careful at "superficial/deep"!
  accuracy.data$MRI.minv[startswith(accuracy.data$MRI.minv,"sup",ignore.case=TRUE)]<-"superficial"
  accuracy.data$MRI.minv[startswith(accuracy.data$MRI.minv,"deep",ignore.case=TRUE)]<-"deep"
  table(accuracy.data$MRI.minv)
  
  accuracy.data<-subset(accuracy.data, histo.minv %in% c("superficial","deep") & MRI.minv %in% c("superficial","deep"))
  table(accuracy.data$histo.minv)
  table(accuracy.data$MRI.minv)
}

print.incorrect.stage<-function(){
  stage.diff <- which(accuracy.data$histo.stage!=accuracy.data$MRI.stage)
  length(stage.diff) # 126 incorrectly staged
  stage.diff
  sapply(stage.diff, function(x) paste0("MRI staged ", accuracy.data$MRI.stage[x], ", whereas tissue examination/histology showed ", accuracy.data$histo.stage[x]))
}

print.incorrect.minv<-function(){
  minv.diff <- which(accuracy.data$histo.minv!=accuracy.data$MRI.minv)
  length(minv.diff) # 73 incorrectly assessed depth of myom. inv.
  sapply(minv.diff, function(x) paste0("invassion depth assessed with MRI as ", accuracy.data$MRI.minv[x], ", whereas tissue examination/histology showed ", accuracy.data$histo.minv[x], " invasion"))
}

# add column from the large dataframe called "data" to another smaller dataframe
add.entire.column<-function(df,raw.name,new.name){
  # extract column index by column name
  cindex<-grep(raw.name,colnames(data))
  # add new column to the dataframe with the name "new.name"
  if(nrow(df)!=0)
  {
    df <- cbind(df, new.name = data[,cindex])
    cnames<-head(colnames(df),-1)
    colnames(df)<-c(cnames,new.name)
  }
  else
  {
    df <- data.frame(new.name = data[,cindex])
    colnames(df)<-c(new.name)
  }
  return(df)
}

add.column<-function(df,raw.name,new.name,rindices){
  cindex<-grep(raw.name,colnames(data))
  if(nrow(df)==0)
  {
    df <-data.frame(data[rindices,cindex])
    colnames(df)<-c(new.name)
  }
  else{
    df <- cbind(df, data[rindices,cindex])
    cnames<-head(colnames(df),-1)
    colnames(df)<-c(cnames,new.name)
  }
  return(df)
}

add.patient.report.data<-function(rind){
  maindf<-data.frame()
  maindf<-add.column(maindf,"PatientID","PatientID",rind)
  maindf
  maindf<-add.column(maindf,"radiomicsCaseID","RadiomicsCaseID",rind)
  maindf
  maindf<-add.column(maindf,"AgeAtDiagnosis","AgeAtDiagnosis",rind)
  maindf
  return(maindf)
}

create.table.rs<-function(ind){
  wrong.risks<-add.patient.report.data(ind)
  wrong.risks<-add.column(wrong.risks,"OpDate","OpDate",ind)
  wrong.risks<-add.column(wrong.risks,"STAGE","HistoStage",ind)
  wrong.risks<-add.column(wrong.risks,"Grade","Grade",ind)
  wrong.risks<-add.column(wrong.risks,"X.1","HistoType",ind)
  wrong.risks<-add.column(wrong.risks,"LVSI","LVSI",ind)
  wrong.risks<-add.column(wrong.risks,"Peritoneal.cyto","PeritonealCytoEval",ind)
  wrong.risks<-add.column(wrong.risks,"Risk.score","RiskScore",ind)
  wrong.risks<-add.column(wrong.risks,"Sup.Deep.2","MyometrialInvasion",ind)
  start<-grep("Y.N.5", colnames(data))
  end<-grep("Y.N.23", colnames(data))
  inv.names<-c("Other involvement: CervicalStroma","Vagina/Parametria","Adnexa","Serosa","Peritoneum/Omentum","R-pelvicLN","L-pelvicLN","Para-aorticLN","R-inguinalLN","L-inguinalLN")
  wrong.risks[inv.names]<-data[ind,seq(start,end,by=5)]
  return(wrong.risks)
}

create.ExcelTable.file<-function(df,file.name)
{
  library(writexl)
  write_xlsx(df,paste0("C:\\Users\\Casi\\Documents\\EndoProj\\",file.name,".xlsx"))
}

clean.risk.scores<-function(){
  Lhisto.findings<-which(data$Risk.score=="low"&(data$X.1=="mucinous"&data$STAGE=="IA"&data$Grade==1 |
                                                   data$X.1=="mucinous"&data$STAGE=="IA"&data$Grade==2 |
                                                   data$X.1=="endometrioid"&data$STAGE=="IA"&data$Grade==1 |
                                                   data$X.1=="endometrioid"&data$STAGE=="IA"&data$Grade==2))
  low.risk<-which(data$Risk.score=="low")
  cat(length(Lhisto.findings),"out of",length(low.risk),"low risk assessments correspond to the histological findings",sep=" ")
  LRdiscrep<-which(data$Risk.score=="low"&!(data$X.1=="mucinous"&data$STAGE=="IA"&data$Grade==1 |
                                              data$X.1=="mucinous"&data$STAGE=="IA"&data$Grade==2 |
                                              data$X.1=="endometrioid"&data$STAGE=="IA"&data$Grade==1 |
                                              data$X.1=="endometrioid"&data$STAGE=="IA"&data$Grade==2))
  #negatia aceastra presupune ca una din fiecare sa fie adevarata
  #the negaties implies that at least one relation from each ANDed line is true
  length(LRdiscrep)
  LRdiscrep
  #LRdiscrep misses one index because it does not know to compare to NA
  #NA!="IA" gives NA
  data$PatientID[setdiff(low.risk,Lhisto.findings)]
  
  all.lrsc<-which(data$Risk.score=="low"&data$Reasons.for.Exclusion %in% c("include","missingMRI"))
  lrsc.histo<-intersect(all.lrsc,which(data$STAGE=="IA"&grepl("sup",data$Sup.Deep.2,ignore.case=TRUE)&!grepl("\\?",data$Sup.Deep.2)&data$X.1=="endometrioid"&data$Grade %in% c(1,2)&(grepl("n",data$LVSI)|data$LVSI=="NO")))
  withoutLVSI<-intersect(all.lrsc,which(data$STAGE=="IA"&grepl("sup",data$Sup.Deep.2,ignore.case=TRUE)&!grepl("\\?",data$Sup.Deep.2)&data$X.1=="endometrioid"&data$Grade %in% c(1,2)))

  length(lrsc.histo)
  length(all.lrsc)
  length(withoutLVSI)
  lrsc.discrep<-setdiff(all.lrsc,lrsc.histo)
  #for lrsc.histo
  table(data$LVSI[lrsc.discrep]) #19 cases non-negative
  table(data$X.1[lrsc.discrep]) #3 mucinous
  table(data$STAGE[lrsc.discrep]) # all IA
  table(data$Grade[lrsc.discrep]) # all grade 1 and 2
  table(data$Sup.Deep.2[lrsc.discrep]) #23 non-sup
  
  
  table(data$LVSI[Lhisto.findings])
  table(data$X.1[Lhisto.findings])
  # 102, 229 - residual disease
  
  
  Hhisto.findings<-which(data$Risk.score=="high"& (data$STAGE=="II" |
                                                     data$STAGE=="IB" & data$Grade==3 & data$X.1=="endometrioid" |
                                                     data$X.1!="endometrioid" & data$X.1!=""))
  high.risk<-which(data$Risk.score=="high")
  #HRdiscrep<-which(data$Risk.score=="high"&data$STAGE!="II")
  cat(length(Hhisto.findings),"out of",length(high.risk),"high risk assessments correspond to the histological findings",sep=" ")
  data$PatientID[setdiff(high.risk,Hhisto.findings)]
  
  
  Ihisto.findings<-which(data$Risk.score=="intermediate"&
                           (data$STAGE=="IA" & data$Grade==3 & data$X.1=="endometrioid" |
                              data$STAGE=="IB" & data$Grade==1 & data$X.1=="endometrioid" |
                              data$STAGE=="IB" & data$Grade==2 & data$X.1=="endometrioid"))
  intermediate.risk<-which(data$Risk.score=="intermediate")
  cat(length(Ihisto.findings),"out of",length(intermediate.risk),"intermediate risk assessments correspond to the histological findings",sep=" ")
  data$PatientID[setdiff(intermediate.risk,Ihisto.findings)]
  
  
  library(sjmisc)
  Ahisto.findings<-which(data$Risk.score=="advanced" & str_contains(data$STAGE,"III") & str_contains(data$STAGE,"IV"))
  Ahisto.findings<-which(data$Risk.score=="advanced" & (str_contains(data$STAGE,"III") | str_contains(data$STAGE,"IV")))
  advanced.risk<-which(data$Risk.score=="advanced")
  cat(length(Ahisto.findings),"out of",length(advanced.risk),"advanced risk assessments correspond to the histological findings",sep=" ")
  data$STAGE[Ahisto.findings]
  
  table(data$STAGE)
  
  rs.ind<-c(setdiff(low.risk,Lhisto.findings),setdiff(high.risk,Hhisto.findings),setdiff(intermediate.risk,Ihisto.findings))
  dff<-create.table.rs(rs.ind)
  dff
  create.ExcelTable.file(dff,"Risk_scores_discrepancies")
  
  wrong.histo.types<-c("no cancer", "uterine metastasis from other cancer", "sarcoma", "stromal", "no histology", "typical HP", "atypical HP", "neuroendocrine", "")
  xhisto.findings<-which(data$Risk.score=='x' & data$X.1 %in% wrong.histo.types)
  missing.risk<-which(data$Risk.score=='x')
  length(xhisto.findings)
  length(missing.risk)
  data$PatientID[setdiff(missing.risk,xhisto.findings)]

  # debug create.table.rs call (add.column function specifically)
  # df1<-data.frame()
  # df1
  # cind<-grep("PatientID",colnames(data))
  # cind
  # nrow(df1)
  # data[rs.ind,cind]
  # if(nrow(df1)!=0){
  #   df1 <<- cbind(df1, "PatientID" = data[rs.ind,cind])
  #   print(df1)
  #   print("if")
  # } else{
  #   df1<<-data.frame("PatientID" = data[rs.ind,cind])
  #   print(df1)
  #   print("else")
  # }
  
  
}

get.MRI.findings<-function(){
  start<-grep("Y.N.4", colnames(data))
  end<-grep("Y.N.22", colnames(data))
  return(seq(start,end,by=5))
}

create.table.study.pop<-function(intind){
  PatientID<-data$PatientID[intind]
  RadiomicsCaseID<-data$radiomicsCaseID[intind]
  DeathDate<-data$Death.Date[intind]
  Reasons.for.Exclusion<-data$Reasons.for.Exclusion[intind]
  intdf<-data.frame(PatientID,RadiomicsCaseID,DeathDate,Reasons.for.Exclusion)
  start<-grep("Y.N.4", colnames(data))
  end<-grep("Y.N.22", colnames(data))
  intdf$"Myometrial invasion"<-data$Depth.of.Myometrial.Invasion[intind]
  cni<-c("Other involvement: Cervical stroma","Vagina/Parametria","Adnexa","Serosa","Peritoneum/Omentum","R-pelvic LN","L-pelvic LN","Para-aortic LN","R-inguinal LN","L-inguinal LN")
  intdf[cni]<-data[intind,seq(start,end,by=5)]
  intdf$"MRI Stage"<-data$X.2[intind]
  return(intdf)
}

get.all.death.info<-function()
{
  noi<-which(data$Death.Date %in% c("no","No","no ") & data$Reasons.for.Exclusion=="include")
  noii<-which(data$Death.Date %in% c("no","No") & data$Reasons.for.Exclusion=="include")
  #slashi<-which(grepl("/",data$Death.Date) & data$Reasons.for.Exclusion=="include")
  #which(grepl("-",data$Death.Date) & data$Reasons.for.Exclusion=="include")
  zeroi<-which(grepl("0",data$Death.Date) & data$Reasons.for.Exclusion=="include")
  #setdiff(slashi,zeroi)
  #159 excluded because of risk score "x" (which she said she checked, so could not be saved)
  all.dd<-c(noii,zeroi)
  #nod it is exactly the same with ndformat
  nod<-which(data$Risk.score %in% c("low","intermediate","high","advanced") & data$Reasons.for.Exclusion %in% c("include","missingMRI") & !grepl("0",data$Death.Date) & data$Death.Date %nin% c("no","No"))
  length(nod)
  him<-setdiff(dataEC416$PatientID,data$PatientID[all.dd])
  me<-setdiff(data$PatientID[all.dd],dataEC416$PatientID)
  setdiff(whatIHave,me)
  setdiff(me,whatIHave)
  setdiff(whatExtHas,him)
  setdiff(him,whatExtHas)
  himid<-which(data$PatientID %in% him)
  ndformat<-which(data$Risk.score %in% c("low","intermediate","high","advanced") & data$Reasons.for.Exclusion %in% c("include","missingMRI") & !grepl("0",data$Death.Date) & data$Death.Date %nin% c("no","No"))
  himnd<-intersect(himid,ndformat)
  setdiff(himnd,nod)
  setdiff(himnd,ndformat)
  setdiff(nod,himnd) #len=62-56=6 to add for missing death info
  rowsToAdd6<-setdiff(ndformat,himnd)
  df6incompl<-create.table.study.pop(rowsToAdd6)
  df6<-add.column(df6incompl,"OpDate","OpDate",rowsToAdd6)
  df6<-add.column(df6,"MRI.date","MRI Date",rowsToAdd6)
  df6<-add.column(df6,"Outcome","Outcome",rowsToAdd6)
  df6<-add.column(df6,"Still.on.FU","Still on FU",rowsToAdd6)
  df6<-add.column(df6,"ImageVoxelSize.mm.","ImageVoxelSize (mm)",rowsToAdd6)
  df6<-add.column(df6,"TumourVolume.mm.3.","TumourVolume (mm^3)",rowsToAdd6)
  df6<-add.column(df6,"N.Voxels","N Voxels",rowsToAdd6)
  df6<-add.column(df6,"used4Radiomics","Used for Radiomics",rowsToAdd6)
  library(writexl)
  write_xlsx(df6,"C:\\Users\\Casi\\Documents\\EndoProj\\Additional_6patients_NoInfoDeathDate.xlsx")
  
  
  noMRIf<-which(data$Risk.score %in% c("low","intermediate","high","advanced") & data$Reasons.for.Exclusion %in% c("missingMRI"))
  himnoMRI<-intersect(himid,noMRIf)
  setdiff(himnoMRI,noMRIf)
  setdiff(noMRIf,himnoMRI)#len=159-89=70 to add for missing MRI data
  rowsToAdd70<-setdiff(noMRIf,himnoMRI)
  df70incompl<-create.table.study.pop(rowsToAdd70)
  df70<-add.column(df70incompl,"OpDate","OpDate",rowsToAdd70)
  df70<-add.column(df70,"MRI.date","MRI Date",rowsToAdd70)
  df70<-add.column(df70,"Outcome","Outcome",rowsToAdd70)
  df70<-add.column(df70,"Still.on.FU","Still on FU",rowsToAdd70)
  df70<-add.column(df70,"ImageVoxelSize.mm.","ImageVoxelSize (mm)",rowsToAdd70)
  df70<-add.column(df70,"TumourVolume.mm.3.","TumourVolume (mm^3)",rowsToAdd70)
  df70<-add.column(df70,"N.Voxels","N Voxels",rowsToAdd70)
  df70<-add.column(df70,"used4Radiomics","Used for Radiomics",rowsToAdd70)
  library(writexl)
  write_xlsx(df70,"C:\\Users\\Casi\\Documents\\EndoProj\\Additional_70patients_missingMRIdata.xlsx")
  
  
  ndformat1<-which(data$Reasons.for.Exclusion %in% c("include","missingMRI") & !grepl("0",data$Death.Date) & data$Death.Date %nin% c("no","No"))
  himnd1<-intersect(himid,ndformat1)
  setdiff(himnd1,ndformat1)
  setdiff(ndformat1,himnd1)#len=62-56=6 to add for missing death info (risk x included)

  noMRIf1<-which(data$Reasons.for.Exclusion %in% c("missingMRI"))
  himnoMRI1<-intersect(himid,noMRIf1)
  setdiff(himnoMRI1,noMRIf1)
  setdiff(noMRIf1,himnoMRI1)#len=162-89=73 to add for missing MRI data  
  setdiff(noMRIf1,noMRIf)
  data$Risk.score[setdiff(noMRIf1,noMRIf)] #excluded bc no risk score
}

all.MRI.recheck<-function(){
  colnames(data30MRI)[1]<-"Radiomics.ID"
  data$PatientID[which(data$radiomicsCaseID %in% data30MRI$Radiomics.ID)]
  data$Reasons.for.Exclusion[which(data$radiomicsCaseID %in% data30MRI$Radiomics.ID)]
  ar.cases<-which(data$radiomicsCaseID %in% data30MRI$Radiomics.ID)
  data$Death.Date[ar.cases]
  data30MRI$DOD
  data$radiomicsCaseID[ar.cases]
  data30MRI$Radiomics.ID
  nointerest<-setdiff(ar.cases,ndformat)#she added some cases I am not interested in
  data$Death.Date[nointerest]#I have the death date for them
  #let's see if it is the same
  dodIhave<-which(data30MRI$Radiomics.ID %in% data$radiomicsCaseID[nointerest])
  big.dataset<-data[order(data$radiomicsCaseID),]
  ar.dataset<-data30MRI[order(data30MRI$Radiomics.ID),]
  nointerest1<-which(big.dataset$radiomicsCaseID %in% ar.dataset$Radiomics.ID[dodIhave])
  Death.Date<-big.dataset$Death.Date[nointerest1]
  DOD<-ar.dataset$DOD[dodIhave]
  radiomicsCaseID<-big.dataset$radiomicsCaseID[nointerest1]
  Radiomics.ID<-ar.dataset$Radiomics.ID[dodIhave]
  compData<-data.frame(Death.Date,DOD,radiomicsCaseID,Radiomics.ID)
  
  
  setdiff(ndformat,ar.cases)
  interesting<-intersect(ndformat,ar.cases)
  data$Death.Date[interesting]#for these I do not have the death date
  #let's see if she gave it to me
  #1 date extra - the radiology findings do not correspond to the MRI findings, should I have them?
  dod.missed<-which(data30MRI$Radiomics.ID %in% data$radiomicsCaseID[interesting])
  data30MRI$DOD[dod.missed]
  
  
  #only numbers
  data$radiomicsCaseID[which(grepl("^[[:digit:]]+$",data$radiomicsCaseID))]
  #the rest of the patients (in this case it is just a special case that the rest contains only letters)
  data$radiomicsCaseID[which(!grepl("^[[:digit:]]+$",data$radiomicsCaseID))]
  #not working: data$PatientID[which(grepl("^[A-Za-z]+$",data$radiomicsCaseID[1]))]#only letters
}

add.categories<-function(){
  library(sjmisc)
  #the constant vector first used in the function clean.risk.scores
  data$"Accepted Histology Type"<-ifelse(data$X.1 %in% wrong.histo.types,"no",data$X.1)
  add.column.ExcelTable("Accepted Histology Type")
  
  other.ca<-c("anal ca","breast cancer","colorectal ca","colon cancer","renal cancer","lung cancer","ovarian cancer")
  data$"Other Malignancies"<-ifelse(grepl(paste0(other.ca,collapse="|"),data$Reasons.for.Exclusion),data$Reasons.for.Exclusion,"none")
  add.column.ExcelTable("Other Malignancies")
  
  MRIfdata<-data[,get.MRI.findings()]
  MRIfdata
  validMRIf<-vector()
  for(i in get.MRI.findings()){
    validMRIf<-c(validMRIf,which(tolower(data[,i]) %in% c("y","n")))
  }
  
  validMRIf.corrected<-vector()
  for(i in get.MRI.findings()){
    validMRIf.corrected<-which(tolower(data[,i]) %in% c("y","n"))
  }
  missing.MRIf.corrected<-setdiff(1:624,validMRIf.corrected)
  
  validMRIf
  missing.MRIf<-setdiff(1:624,validMRIf)
  missing.MRIf
  MRI.minv <- ifelse(data$"Sup.Deep"!="x" & data$"Sup.Deep"!="", data$"Sup.Deep", data$"Sup.Deep.1")
  data$"Depth of Myometrial Invasion"<-MRI.minv
  add.column.ExcelTable("Depth of Myometrial Invasion")
  
  valid.MRIminv<-which(grepl("sup",MRI.minv,ignore.case = TRUE)&!grepl("deep",MRI.minv,ignore.case=TRUE)|
                       grepl("deep",MRI.minv,ignore.case = TRUE)&!grepl("sup",MRI.minv,ignore.case=TRUE))
  MRI.minv
  valid.MRIminv
  missing.MRIminv<-setdiff(1:624,valid.MRIminv)
  missing.MRIminv
  table(MRI.minv)
  table(data$X.2)
  missing.MRIst<-which(data$X.2=="")
  missing.MRIst
  data$"Availability of MRI"<-"yes"
  data$`Availability of MRI`
  data$`Depth of Myometrial Invasion`
  data$"Availability of MRI"[missing.MRIf]<-"missing y/n in MRI findings"
  
  free.indices<-which(data$`Availability of MRI`=="yes")
  mindices<-intersect(free.indices,missing.MRIminv)
  extra.mindices<-setdiff(missing.MRIminv,mindices)
  data$"Availability of MRI"[mindices]<-"missing MRI depth of invasion ONLY"
  data$"Availability of MRI"[extra.mindices]<-"missing depth of invasion and other MRI findings"
  
  stindices<-intersect(which(data$`Availability of MRI`=="yes"),missing.MRIst)
  ynind<-intersect(which(grepl("y/n",data$`Availability of MRI`)),missing.MRIst)
  minv.only.ind<-intersect(which(grepl("ONLY",data$`Availability of MRI`)),missing.MRIst)
  morefind<-intersect(which(grepl("other",data$`Availability of MRI`)),missing.MRIst)
  data$"Availability of MRI"[stindices]<-"missing MRI stage ONLY"
  data$"Availability of MRI"[ynind]<-"missing MRI stage but also y/n in MRI findings"
  data$"Availability of MRI"[minv.only.ind]<-"missing MRI stage and depth of invasion"
  data$"Availability of MRI"[morefind]<-"missing MRI stage, depth of invasion, and other MRI findings"
  add.column.ExcelTable("Availability of MRI")
  # data$"Availability of MRI"[intersect(missing.MRIf,missing.MRIminv)]<-"missing depth of invasion and other MRI findings"
  # data$"Availability of MRI"[intersect(missing.MRIf,missing.MRIst)]<-"missing MRI stage and depth of invasion"
  # data$"Availability of MRI"[intersect(missing.MRIminv,missing.MRIst)]<-"missing MRI stage and depth of invasion"
  tail(data)
  
  # risk scores with regard to the patient outcomes using a Kaplan-Meier curve
  incl.ind<-intersect(which(grepl("missing",data$Availability.of.MRI)),which(data$Reasons.for.Exclusion=="include"))
  data$Reasons.for.Exclusion[incl.ind]<-"missingMRI"
  data$Availability.of.MRI[which(data$Reasons.for.Exclusion=="missingMRI")]
  
  table(data$Reasons.for.Exclusion)
  which(data$Availability.of.MRI!="yes")
  
  colnames(data)[130]<-"Reasons for Exclusion"
  
  #adding missingMRI reason to already excluded patients
  `%nin%` = Negate(`%in%`)
  alr.excl.ind<-which(grepl("missing",data$Availability.of.MRI) & !grepl("missing",data$Reasons.for.Exclusion))
                          #which(data$Reasons.for.Exclusion %nin% c("missingMRI","include")))
  data$Reasons.for.Exclusion[alr.excl.ind]
  #data$Reasons.for.Exclusion[which(data$Reasons.for.Exclusion %nin% c("missingMRI","include"))]
  data$Reasons.for.Exclusion[alr.excl.ind]<-add.reason("@ the beginning",alr.excl.ind,"missingMRI")

  
  #colnames(dataEC416)[1]<-"PatientID"
  dataEC416$PatientID
  inclInStudy<-which(data$Reasons.for.Exclusion=="include" &
  (grepl("0",data$Death.Date) | data$Death.Date %in% c("no","No","no ")))#, "ND")))
  #0 instead of -
  #which(data$Reasons.for.Exclusion=="include" &
          #(grepl("/",data$Death.Date) | data$Death.Date %in% c("no","No","no ")))
  whatExtHas<-setdiff(dataEC416$PatientID,data$PatientID[inclInStudy])
  whatIHave<-setdiff(data$PatientID[inclInStudy],dataEC416$PatientID)
  intsct<-intersect(data$PatientID[inclInStudy],dataEC416$PatientID)
  pind<-which(data$PatientID %in% intsct)#get indexes from ids
  commind<-which(inclInStudy %in% pind)#det the position the common patients are in my list
  os.time[commind]
  commind.ext<-which(dataEC416$PatientID %in% intsct)
  dataEC416$"OpDate2Aug3rd2020.day."[commind.ext]
  difftime("2020-08-03","2012-02-08",tz="GMT", units = "days")
  difftime("2017-03-29","2012-02-09",tz="GMT", units = "days")
  difftime("2020-08-03","2012-02-09",tz="GMT", units = "days")
  #he considered 3 august 2020(=the censorship date) even for the patients 
  #who died before the end of the study
  
  extind<-which(data$PatientID %in% whatExtHas)
  data$Death.Date[extind]#they are not ordered in ext table
  rformat<-which(grepl("-",data$Death.Date) | data$Death.Date %in% c("no","No","no "))
  # patients from the ext table with the 
  #right format for Death Date, but with missing MRI
  table(data$Reasons.for.Exclusion[intersect(extind,rformat)])
  # for the remaining patients in external table, 
  #there is missing info on Death Date 
  table(data$Death.Date[setdiff(extind,rformat)])
  #whynot<-which(data$Death.Date[extind] %in% c("no","No")|grepl("-",data$Death.Date[extind]))
  #data$Death.Date[whynot]#not good because index is of data$Death.Date[extind]
  setdiff(whatExtHas,data$PatientID[extind])#by comparing both you check you got the right indices from EndoRad nbs
  
  PatientID<-data$PatientID[extind]
  RadiomicsCaseID<-data$radiomicsCaseID[extind]
  DeathDate<-data$Death.Date[extind]
  Reasons.for.Exclusion<-data$Reasons.for.Exclusion[extind]
  extdf<-data.frame(PatientID,RadiomicsCaseID,DeathDate,Reasons.for.Exclusion)
  extdf$Error<-ifelse(!grepl("/",extdf$DeathDate)&extdf$DeathDate %nin% c("no","No","no "),"no information about death","missing MRI data")
  #be careful about dying/no
  extdf$Error<-ifelse(!grepl("/",extdf$DeathDate)&extdf$DeathDate %nin% c("no","No","no ")&extdf$Reasons.for.Exclusion=="missingMRI","both: no information about death & missing MRI data",extdf$Error)
  start<-grep("Y.N.4", colnames(data))
  end<-grep("Y.N.22", colnames(data))
  extdf$"Myometrial invasion"<-data$Depth.of.Myometrial.Invasion[extind]
  cns<-c("Other involvement: Cervical stroma","Vagina/Parametria","Adnexa","Serosa","Peritoneum/Omentum","R-pelvic LN","L-pelvic LN","Para-aortic LN","R-inguinal LN","L-inguinal LN")
  length(cns)
  extdf[cns]<-data[extind,seq(start,end,by=5)]
  extdf$"MRI Stage"<-data$X.2[extind]
  extdf
  library(writexl)
  write_xlsx(extdf,"C:\\Users\\Casi\\Documents\\EndoProj\\Extra_Patients.xlsx")
  
  intind<-which(data$PatientID %in% whatIHave)
  table(data$Reasons.for.Exclusion[intind])
  table(data$Death.Date[intind])
  PatientID<-data$PatientID[intind]
  RadiomicsCaseID<-data$radiomicsCaseID[intind]
  DeathDate<-data$Death.Date[intind]
  Reasons.for.Exclusion<-data$Reasons.for.Exclusion[intind]
  intdf<-data.frame(PatientID,RadiomicsCaseID,DeathDate,Reasons.for.Exclusion)
  start<-grep("Y.N.4", colnames(data))
  end<-grep("Y.N.22", colnames(data))
  intdf$"Myometrial invasion"<-data$Depth.of.Myometrial.Invasion[intind]
  cni<-c("Other involvement: Cervical stroma","Vagina/Parametria","Adnexa","Serosa","Peritoneum/Omentum","R-pelvic LN","L-pelvic LN","Para-aortic LN","R-inguinal LN","L-inguinal LN")
  length(cni)
  intdf[cni]<-data[intind,seq(start,end,by=5)]
  intdf$"MRI Stage"<-data$X.2[intind]
  intdf
  library(writexl)
  write_xlsx(intdf,"C:\\Users\\Casi\\Documents\\EndoProj\\Extra_Patients.xlsx")
  
  #incomplete!!
  nd<-which(data$Reasons.for.Exclusion=="include"&!grepl("/",data$Death.Date) & data$Death.Date %nin% c("no","No","no "))
  #all included cases without death info except for "dying/no"
  dataf<-data.frame(data$PatientID[nd],data$AgeAtDiagnosis[nd],data$radiomicsCaseID[nd],data$MRI.date[nd],data$Date.MDT[nd],data$OpDate[nd],data$Last.seen[nd],data$Still.on.FU[nd],data$Death.Date[nd])
  colnames(dataf)<-c("Patient ID","Age at Diagnosis","Radiomics Case ID","MRI Date","MDT Date","Surgery Date","Last Seen","Still on FU","Death Date")
  dataf$"Any Reasons for Exclusion?"<-data$Reasons.for.Exclusion[nd]
  #dataf$"Days of Survival (from op.date to death.date)"<-
  #dataf$"Days Sx to RIP"<-
  #dataf$"Event Status"<-ifelse(os.event[nd]==1,"dead","alive")
  library(writexl)
  write_xlsx(dataf,"C:\\Users\\Casi\\Documents\\EndoProj\\NoInfo_on_DeathDate.xlsx")
  
  
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
  
  
  condition<-which(grepl("recurrence",dataS2$Date.op,ignore.case=TRUE) | 
   grepl("recurrence",dataS2$Notes,ignore.case=TRUE) | 
   grepl("recurrence",dataS2$X,ignore.case=TRUE))
  ids<-dataS2$PatientID[condition]
  rn<-which(data$PatientID %in% id)
  #rn<-intersect(583:627,which(...))
  data$Reasons.for.Exclusion[rn]<- add.reason("@ the end",rn, "recurrence")
  
  conditionInfo<-which(grepl("info",dataS2$X,ignore.case=TRUE) | 
                      grepl("No note on cerner",dataS2$Notes,ignore.case=TRUE) |
                      grepl("Nothing",dataS2$Notes,ignore.case=TRUE) |
                      dataS2$Date.op==""&dataS2$X=="")#!grepl("recurrent",dataS2$Notes,ignore.case=TRUE)&!grepl("benign",dataS2$Notes,ignore.case=TRUE)))
  idsInfo<-dataS2$PatientID[conditionInfo]
  rnInfo<-which(data$PatientID %in% idsInfo)
  data$Reasons.for.Exclusion[rnInfo]<-add.reason("@ the end",rnInfo, "insufficient surgical information")
  
  conditionAdv<-which(grepl("advanced",dataS2$Date.op,ignore.case=TRUE) |
                      grepl("advanced",dataS2$X,ignore.case=TRUE))
  idsAdv<-dataS2$PatientID[conditionAdv]
  rnAdv<-which(data$PatientID %in% idsAdv)
  data$Reasons.for.Exclusion[rnAdv]<-add.reason("@ the end",rnAdv, "no surgery (advanced disease)")
  
  conditionElsew<-which(grepl("elsewhere",dataS2$Date.op,ignore.case=TRUE) |
                        grepl("elsewhere",dataS2$X,ignore.case=TRUE) |
                        grepl("operated on a diff hospital",dataS2$X,ignore.case=TRUE))
  idsElsew<-dataS2$PatientID[conditionElsew]
  rnElsew<-which(data$PatientID %in% idsElsew)
  data$Reasons.for.Exclusion[rnElsew]<-add.reason("@ the end",rnElsew, "surgery elsewhere")
  #which(grepl("surgery elsewhere",data$Reasons.for.Exclusion))
  
  particase1<-which(data$PatientID=="EndoRad_616")
  data$Reasons.for.Exclusion[particase1]<-paste0(data$Reasons.for.Exclusion[particase1]," & no surgery (treated conservatively elsewhere)")
  
  particase2<-which(data$PatientID=="EndoRad_622")
  data$Reasons.for.Exclusion[particase2]<-paste0(data$Reasons.for.Exclusion[particase2]," & no cancer")
  
  table(dataS2$Palliative)
  table(dataS2$Palliative.)#Palliative? this one is accurate!
  conditionPalliv<-which(grepl("palliative",dataS2$X,ignore.case=TRUE)|grepl("palliative",dataS2$Notes,ignore.case=TRUE)
                         | dataS2$Palliative.=="y")
  
  conditionCons<-which(!grepl("young",dataS2$Date.op,ignore.case=TRUE)&!grepl("Dubai",dataS2$Notes)&(grepl("conservative",dataS2$X,ignore.case=TRUE) |
                       grepl("conservative",dataS2$Date.op,ignore.case=TRUE)))
  idsCons<-dataS2$PatientID[conditionCons]
  rnCons<-which(data$PatientID %in% idsCons)
  data$Reasons.for.Exclusion[rnCons]<-add.reason("@ the end",rnCons, "no surgery (treated conservatively)")
  
  conditionYoung<-which(grepl("young",dataS2$Date.op,ignore.case=TRUE))
  idsYoung<-dataS2$PatientID[conditionYoung]
  rnYoung<-which(data$PatientID %in% idsYoung)
  data$Reasons.for.Exclusion[rnYoung]<-add.reason("@ the end",rnYoung, "too young for surgery (treated conservatively)")
  
  #remainingIds<-setdiff(1:49,allcond)
  all_conditionProvera<-which(grepl("provera",dataS2$Notes,ignore.case=TRUE) |
                          grepl("progesterone therapy",dataS2$Notes))
  #conditionProvera<-intersect(remainingIds,all_conditionProvera)
  conditionProvera<-setdiff(all_conditionProvera,c(25,33,41))
  idsProvera<-dataS2$PatientID[conditionProvera]
  rnProvera<-which(data$PatientID %in% idsProvera)
  data$Reasons.for.Exclusion[rnProvera]<-add.reason("@ the end",rnProvera, "on Provera")
  
  conditionDecl<-which(grepl("decline",dataS2$S.P.D.) | 
        grepl("decline",dataS2$Notes,ignore.case=TRUE))
  
  all_conditionFit<-which(grepl("fit",dataS2$Date.op,ignore.case=TRUE) |
                      grepl("fit",dataS2$X,ignore.case=TRUE))
  # 10 is not fit, but I wrote it manually
  conditionFit<-setdiff(all_conditionFit,c(8,10,15))
  idsFit<-dataS2$PatientID[conditionFit]
  rnFit<-which(data$PatientID %in% idsFit)
  data$Reasons.for.Exclusion[rnFit]<-add.reason("@ the end",rnFit, "not fit for surgery")
  
  
  allcond<-c(condition,conditionInfo,conditionAdv,conditionElsew,41,47,28,conditionCons,conditionYoung,15,all_conditionProvera,conditionDecl,44,conditionFit)
  setdiff(1:49,allcond)
  
  
  add.column.ExcelTable("Reasons.for.Exclusion")
  data$Reasons.for.Exclusion
}

add.reason<-function(place, pos, word)
{
  if(place=="@ the beginning")
    #paste0("missingMRI & ",data$Reasons.for.Exclusion[alr.excl.ind])
    sth<-sapply(data$Reasons.for.Exclusion[pos],function(x){x<-paste0(word," & ",x)},USE.NAMES = FALSE)
  else if(place=="@ the end")
    sth<-sapply(data$Reasons.for.Exclusion[pos],function(x){x<-paste0(x," & ",word)},USE.NAMES = FALSE)
  return(sth)
 
  #alternatives: lapply and 
  #data$Reasons.for.Exclusion[row.nb]<-mapply(paste0,data$Reasons.for.Exclusion[row.nb]," & recurrence",USE.NAMES = FALSE)
  #mapply(paste0,"ana","are",USE.NAMES = FALSE)
}

# time between diagnosis and surgical treatment - demographic data
data$MRI.date[allInclStudy]
intersect(which(data$MRI.date==""),allInclStudy)
#EndoRad_430 all MRI data but no MRI date!! 
opd<-sub(" .*", "", data$OpDate[allInclStudy[-213]])
#both formatting data in that order
#format(as.Date(data$MRI.date[allInclStudy[-213]]), "%Y-%m-%d")
#strftime(data$MRI.date[allInclStudy[-213]], "%Y-%m-%d")
MRId<-strftime(strptime(data$MRI.date[allInclStudy[-213]],"%d/%m/%Y"),"%Y-%m-%d")
difftime(opd,MRId, tz="GMT", units = "days")

#part of data cleaning
data$Grade[which(data$X.1 %in% c("clear cell","serous","carcinosarcoma"))]
table(data$Grade[which(data$X.1=="clear cell")])
table(data$Grade[which(data$X.1=="serous")])
id2correct<-which(data$X.1=="serous" & data$Grade=="2")
data$Reasons.for.Exclusion[id2correct]
corrected.data$Reasons.for.Exclusion[id2correct]
#modification changed manually in the table
corrected.data$Grade[id2correct]<-"3"
table(corrected.data$Grade[which(corrected.data$X.1=="serous")])
table(data$Grade[which(data$X.1=="carcinosarcoma")])
table(data$Grade[which(data$X.1=="undifferentiated")])
table(data$Grade[which(data$X.1=="mixed high grade")])

corrected.data$Reasons.for.Exclusion[which(corrected.data$radiomicsCaseID=="396")]
corrected.data$Death.Date[which(corrected.data$radiomicsCaseID=="396")]

#modifications to change in the table
noInvasion<-which(corrected.data$Availability.of.MRI=="missing MRI depth of invasion ONLY" & corrected.data$Depth.of.Myometrial.Invasion=="none")
corrected.data$X.2[intersect(noInvasion,which(corrected.data$Reasons.for.Exclusion=="include"))]
corrected.data$Availability.of.MRI[noInvasion]<-"yes"
corrected.data$Reasons.for.Exclusion[noInvasion]
missingMRIOnly<-intersect(noInvasion, which(corrected.data$Reasons.for.Exclusion=="missingMRI"))
othereasonsAlso<-intersect(noInvasion, grep("breast cancer", corrected.data$Reasons.for.Exclusion))
corrected.data$Reasons.for.Exclusion[missingMRIOnly]<-"include"
corrected.data$Reasons.for.Exclusion[othereasonsAlso]<-"breast cancer"

#modifications to change in the table
table(corrected.data$STAGE)
table(corrected.data$X.2)
corrected.data$STAGE[which(corrected.data$STAGE=="IVB")]<-"IV"
corrected.data$X.2[which(corrected.data$X.2=="IVB")]<-"IV"

#modification to change in the table
corrected.data$Reasons.for.Exclusion[which(corrected.data$radiomicsCaseID=="360")]<-"duplicate"

#AR table
corrected.data$X.2[179]
corrected.data$Custom.worklist[179]

#under surveillance for breast cancer
excluded<-which(corrected.data$radiomicsCaseID=="445")
corrected.data$Reasons.for.Exclusion[excluded]
corrected.data$Availability.of.MRI[excluded]

corrected.data$Reasons.for.Exclusion[which(corrected.data$radiomicsCaseID=="360")]
corrected.data$Availability.of.MRI[which(corrected.data$radiomicsCaseID=="360")]

corrected.data$Reasons.for.Exclusion[which(corrected.data$radiomicsCaseID=="354")]
corrected.data$Availability.of.MRI[which(corrected.data$radiomicsCaseID=="354")]

corrected.data$Reasons.for.Exclusion[which(corrected.data$radiomicsCaseID=="394")]
corrected.data$Availability.of.MRI[which(corrected.data$radiomicsCaseID=="394")]
corrected.data$X.2[which(corrected.data$radiomicsCaseID=="394")]

#not related to NB and AR tables
#dataframe modification to change in the table as well
id2excl<-which(corrected.data$Reasons.for.Exclusion=="include"&corrected.data$Risk.score=="x")
corrected.data$Reasons.for.Exclusion[id2excl]
corrected.data$Risk.score[id2excl]
corrected.data$Reasons.for.Exclusion[id2excl]<-"no risk score assessment"

#NB table - all excluded
corrected.data$Reasons.for.Exclusion[which(corrected.data$radiomicsCaseID=="233")]
corrected.data$Reasons.for.Exclusion[which(corrected.data$radiomicsCaseID=="397")]
corrected.data$Reasons.for.Exclusion[which(corrected.data$radiomicsCaseID=="446")]

corrected.data$Availability.of.MRI[which(corrected.data$radiomicsCaseID=="233")]
corrected.data$Availability.of.MRI[which(corrected.data$radiomicsCaseID=="397")]
corrected.data$Availability.of.MRI[which(corrected.data$radiomicsCaseID=="446")]<-"missing MRI stage but also y/n in MRI findings"

#the following should be included (from NB)
#septate uterus
corrected.data$Reasons.for.Exclusion[which(corrected.data$radiomicsCaseID=="365")]
corrected.data$Reasons.for.Exclusion[which(corrected.data$radiomicsCaseID=="478")]
corrected.data$Reasons.for.Exclusion[which(corrected.data$radiomicsCaseID=="462")]
corrected.data$Availability.of.MRI[which(corrected.data$radiomicsCaseID=="462")]
corrected.data$Availability.of.MRI[which(corrected.data$radiomicsCaseID=="365")]
corrected.data$Availability.of.MRI[which(corrected.data$radiomicsCaseID=="478")]

#ids4 contains all NB except for the ones EXCLUDED
for(i in cids)
  print(table(corrected.data[ids4,i]))
corrected.data[ids4,cids]
all2excl<-vector()
for(i in cids){
  ids2excl<-intersect(ids4,which(grepl("no",corrected.data[,i],ignore.case = TRUE)&nchar(corrected.data[,i])>4))
  print(ids2excl)
  if(length(ids2excl)!=0){
    print("adding indices...")
    all2excl<-c(all2excl,ids2excl)
  }
}
#a pair of indices appears twice because missing info is checked per/by column 
#and 2 patients miss MRI information on two columns 

id2modify<-unique(all2excl)
corrected.data$Reasons.for.Exclusion[id2modify]
corrected.data$Availability.of.MRI[id2modify]
#dataframe modification to change in the table as well
corrected.data$Availability.of.MRI[475]<-"missing MRI depth of invasion ONLY"
corrected.data$Availability.of.MRI[457]<-"missing y/n in MRI findings"
corrected.data$Availability.of.MRI[463]<-"missing y/n in MRI findings"
corrected.data$Availability.of.MRI[473]<-"missing y/n in MRI findings"
corrected.data$Availability.of.MRI[508]<-"missing y/n in MRI findings"

ids2incl<-setdiff(ids4,id2modify)
corrected.data$radiomicsCaseID[ids2incl]
corrected.data$Availability.of.MRI[ids2incl]<-"yes"
corrected.data$Reasons.for.Exclusion[ids2incl]<-"include"

missingMRIf.temp<-vector()
for(i in get.MRI.findings()){
  missingMRIf.temp<-c(missingMRIf.temp,which(tolower(corrected.data[,i]) %nin% c("y","n")))
}
missing.MRIf.toBcorrected<-unique(missingMRIf.temp)
validMRIf.toBcorrected<-setdiff(1:624,missing.MRIf.toBcorrected)

#before assigning a new val to shouldNotHaveMRIf checking if I wrongly filled some fields: 
#compare the old shouldNotHaveMRIf with the new one, setdiff(unique(validMRIf),validMRIf.toBcorrected)

#did I put some additional MRIs, which should not be there?
setdiff(shouldNotHaveMRIf,setdiff(unique(validMRIf),validMRIf.toBcorrected))

#if not, I should just complete the remaining:
remainingShouldNotHaveMRIf<-setdiff(setdiff(unique(validMRIf),validMRIf.toBcorrected),shouldNotHaveMRIf)

#ultimate check: the number of additional indices + the number of indices of the old shouldNotHaveMRIf
#should give the number of indices of the new shouldNotHaveMRIf 
length(remainingShouldNotHaveMRIf)+length(shouldNotHaveMRIf)==length(setdiff(unique(validMRIf),validMRIf.toBcorrected))

#now that we cleared that, we can safely make the assessment:
#change the value of the old to the new one (update shouldNotHaveMRIf)

#keep the old value just in case
shouldNotHaveMRIf1<-shouldNotHaveMRIf
shouldNotHaveMRIf<-setdiff(unique(validMRIf),validMRIf.toBcorrected)
#initially(i.e. when Avail.of.MRI col has been created) wrongly assessed as valid MRI findings data
table(corrected.data$Availability.of.MRI[shouldNotHaveMRIf])
corrected.data[shouldNotHaveMRIf,cids]
corrected.data[remainingShouldNotHaveMRIf,cids]
#for old shouldNotHaveMRIf
corrected.data$Availability.of.MRI[intersect(shouldNotHaveMRIf,which(corrected.data$Availability.of.MRI=="missing MRI stage and depth of invasion"))]<-"missing MRI stage, depth of invasion, and other MRI findings"
corrected.data$Availability.of.MRI[intersect(shouldNotHaveMRIf,which(corrected.data$Availability.of.MRI=="missing MRI stage ONLY"))]<-"missing MRI stage but also y/n in MRI findings"
corrected.data$Availability.of.MRI[intersect(shouldNotHaveMRIf,which(corrected.data$Availability.of.MRI=="yes"))]<-"missing y/n in MRI findings"
corrected.data$Reasons.for.Exclusion[intersect(shouldNotHaveMRIf,which(corrected.data$Availability.of.MRI=="missing y/n in MRI findings"))]<-"missingMRI"
#for the remaining values to change from the new shouldNotHaveMRIf
corrected.data$Availability.of.MRI[intersect(remainingShouldNotHaveMRIf,which(corrected.data$Availability.of.MRI=="yes"))]<-"missing y/n in MRI findings"
corrected.data$Availability.of.MRI[intersect(remainingShouldNotHaveMRIf,which(corrected.data$Availability.of.MRI=="missing MRI stage and depth of invasion"))]<-"missing MRI stage, depth of invasion, and other MRI findings"
corrected.data$Reasons.for.Exclusion[intersect(remainingShouldNotHaveMRIf,which(corrected.data$Availability.of.MRI=="missing y/n in MRI findings"))]<-"missingMRI"

shouldHaveMRIf<-setdiff(missing.MRIf,missing.MRIf.toBcorrected)
#are there patients which have MRI finds and we previously said they have not
corrected.data[shouldHaveMRIf,cids]
table(corrected.data$Availability.of.MRI[shouldHaveMRIf])
#corrected.data$Availability.of.MRI[intersect(shouldHaveMRIf,which(corrected.data$Availability.of.MRI=="missing MRI stage ONLY"))]<-"missing MRI stage but also y/n in MRI findings"
#corrected.data$Availability.of.MRI[intersect(shouldHaveMRIf,which(corrected.data$Availability.of.MRI=="missing MRI stage and depth of invasion"))]<-"missing MRI stage, depth of invasion, and other MRI findings"
#corrected.data$Availability.of.MRI[intersect(shouldHaveMRIf,which(corrected.data$Availability.of.MRI=="yes"))]<-"missing y/n in MRI findings"
corrected.data$Availability.of.MRI[intersect(shouldHaveMRIf,which(corrected.data$Availability.of.MRI=="missing MRI stage, depth of invasion, and other MRI findings"))]<-"missing MRI stage and depth of invasion"
corrected.data$Reasons.for.Exclusion[intersect(shouldHaveMRIf,which(corrected.data$Availability.of.MRI=="missing MRI stage and depth of invasion"))]


library(xlsx)
workbook_vmar15<-loadWorkbook(file = "EC_data_11mar2021_NBlatest.xlsx")
allsheets<-getSheets(workbook_vmar15)
#Ids4<-ids4+3
allrows<-getRows(allsheets$AllCases,rowIndex = 4:627)
ccfin<-getCells(allrows,colIndex = 1:135)
library(stringr)
for(i in 1:length(ccfin)){
  nmi<-names(ccfin)[i]
  vs<-str_extract_all(nmi, "[[:digit:]]+")
  rnb<-as.numeric(vs[[1]][1])-3
  cnb<-as.numeric(vs[[1]][2])
  setCellValue(ccfin[[nmi]],corrected.data[rnb,cnb])
}
saveWorkbook(workbook_vmar15,file = "EC_data_15mar2021.xlsx")

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

# a leap year has 366 days
#NOT USED
years<-seq(2013,2020)
leap.years=vector()
common.years=vector()
for(i in years){
  if(i%%4==0 && i%%100 || i%%400==0)
    leap.years<-c(leap.years,i)
  else
    common.years<-c(common.years,i)
}
#NOT USED
}
#365*6+2*366-the number of days between 2013 and 2020 inclusive
dsurv<-os.time1[-213]

print(n.year.surv(3))
print(n.year.surv(1))
print(n.year.surv(5))
print(n.year.surv(10))
which(dy-opy>1)
n.year.surv<-function(yn,allincl){
  dds<-sub(" .*", "", corrected.data$Death.Date[allincl])
  opd<-sub(" .*", "", corrected.data$OpDate[allincl])
  dd<-ifelse(dds=="no" | dds=="No", "2020-08-03", dds)
  opd[which(dds=="no" | dds=="No")]<-"2007-05-02"
  opy<-as.numeric(substring(opd,1,4))
  dy<-as.numeric(substring(dd,1,4))
  potentials<-which(dy-opy>=yn)
  #allincl<-allInclStudy[-213]
  potential<-which(dy-opy==yn)
  # intersect(which(dy-opy>yn) %in% opd)
  # opd[which(dy-opy>yn)]
  
  #CAN BE IGNORED
  # all.opdates<-sub(" .*", "", data$OpDate)
  # pot.opdates<-opd[potentials]
  # ##which(pot.opdates %in% all.opdates) #toate sunt incluse, bineinteles
  # which(all.opdates %in% pot.opdates)
  # ##intersect(which(all.opdates %in% pot.opdates),allincl)
  # all.ddates<-sub(" .*", "", data$Death.Date)
  # pot.ddates<-dd[potentials]
  # which(all.ddates %in% pot.ddates)
  #CAN BE IGNORED
  
  pid<-corrected.data$PatientID[allincl]
  pot.pid<-pid[potentials]
  allpot.data<-which(corrected.data$PatientID %in% pot.pid)
  
  
  opmon<-as.numeric(substring(opd[potential],6,7))
  dmon<-as.numeric(substring(dd[potential],6,7))
  opday<-as.numeric(substring(opd[potential],9,10))
  dday<-as.numeric(substring(dd[potential],9,10))
  not.nysurvived<-which(dmon<opmon | (dmon==opmon & dday<opday))
  #nleapy<-length(intersect(seq(opy,dy),leap.years))
  #if(dsurv>=nleapy*366+365*(yn-nleapy)
  howManyAlive<-length(potentials)-length(not.nysurvived)
  whoIsAlive<-setdiff(allpot.data,which(corrected.data$PatientID %in% pot.pid[not.nysurvived]))
  nysurvival<-(howManyAlive/length(allincl))*100
  print(whoIsAlive)
  print(nysurvival)
  return(nysurvival)
}

n.year.surv<-function(yn,allincl){
  dds<-sub(" .*", "", corrected.data$Death.Date[allincl])
  opd<-sub(" .*", "", corrected.data$OpDate[allincl])
  dd<-ifelse(dds=="no" | dds=="No", "2020-08-03", dds)
  opy<-as.numeric(substring(opd,1,4))
  dy<-as.numeric(substring(dd,1,4))
  potentials<-which(dy-opy>=yn)
  #allincl<-allInclStudy[-213]
  potential<-which(dy-opy==yn)
  # intersect(which(dy-opy>yn) %in% opd)
  # opd[which(dy-opy>yn)]
  
  pid<-corrected.data$PatientID[allincl]
  pot.pid<-pid[potentials]
  allpot.data<-which(corrected.data$PatientID %in% pot.pid)
  
  
  opmon<-as.numeric(substring(opd[potential],6,7))
  dmon<-as.numeric(substring(dd[potential],6,7))
  opday<-as.numeric(substring(opd[potential],9,10))
  dday<-as.numeric(substring(dd[potential],9,10))
  not.nysurvived<-which(dmon<opmon | (dmon==opmon & dday<opday))
  #nleapy<-length(intersect(seq(opy,dy),leap.years))
  #if(dsurv>=nleapy*366+365*(yn-nleapy)
  howManyAlive<-length(potentials)-length(not.nysurvived)
  whoIsAlive<-setdiff(allpot.data,which(corrected.data$PatientID %in% pot.pid[not.nysurvived]))
  nysurvival<-(howManyAlive/length(allincl))*100
  print(whoIsAlive)
  print(nysurvival)
  return(nysurvival)
}
         
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

analysis<-function(allInclStudy){
  i<-which(risks=="x")#"no ", risk.score=x
  inclInStudy1<-which(data$Reasons.for.Exclusion=="include" &
                       (grepl("0",data$Death.Date) | data$Death.Date %in% c("no","No")))
  allInclStudy<-intersect(inclInStudy1,which(!is.na(data$TumourVolume.mm.3.)))
  allInclStudy1<-allInclStudy
  
  x<-data$TumourVolume.mm.3.[allInclStudy]
  ddsub1<-sub(" .*", "", data$Death.Date[allInclStudy])
  odate1<-sub(" .*", "", data$OpDate[allInclStudy])
  ddate1<-ifelse(ddsub1=="no" | ddsub1=="No", "2020-08-03", ddsub1)
  # strptime(odate1, format = "%Y-%m-%d")
  # strptime(ddate1, format = "%d/%m/%y")
  os.time1<-difftime(ddate1, odate1, tz="GMT", units = "days")
  #format="%dd%/%mm%/%yyyy%"
  os.event1<-as.numeric(ddate1!="2020-08-03")#is.deceased
  #see if there are deaths in the same year as the year of study
  dead<-sum(os.event1==1)
  alive<-sum(os.event1==0)
  total<-length(allInclStudy)#or dead+alive
  dead.pct<-(dead*100)/total
  alive.pct<-(alive*100)/total
  y<-os.time1
  plot(x/100, y, main = "Simple Scatter of Survival (Days) by Tumour Volume (mm^3)",
       xlab = "Tumour Volume (mm^3)", ylab = "Days of Survival",
       pch = 19, frame = FALSE)
}
#add.new.columns
allInclStudy<-intersect(allstudypop,which(!is.na(corrected.data$TumourVolume.mm.3.)))
length(allInclStudy)
allInclStudy2<-allInclStudy

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

agep<-corrected.data$AgeAtDiagnosis[allInclStudy]
library(survival)
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


library(survival)
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


library(survminer)
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

plot(survfit(ec.surv1 ~ stage.MRI1), col=c("black", "red"), fun="cloglog")
plot(survfit(ec.surv1 ~ grade.histo1), col=c("black", "red"), fun="cloglog")
plot(survfit(ec.surv1 ~ strata(grade.histo1)), col=c("black", "red"), fun="cloglog")
#stage.MRI1+age1+tumour.vol1+strata(grade.histo1)
library(ggplot2)
library(survminer)
plot(tumour.vol1,ec.surv1[,1])
plot(age1,ec.surv1[,1])
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
# par(mfrow=c(1,1))
# plot(cox.zph(coxph.model5)[4])#stratify for grade histo!
# abline(h=0,col="red")
# 
# par(mfrow=c(1,1))
# plot(cox.zph(coxph.model5)[5])#stratify for stage histo!
# abline(h=0,col="red")
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

plotIDa.predstrength.survCurvComp<-function(subpop){
  ggls<-list("with MRI stage + CRS + age + tumour.groupingII"=survfit(model11all,subset=subpop),"with CRS + age + tumour.groupingII"=survfit(model11alla,subset=subpop),"with MRI stage + age + tumour.groupingII"=survfit(model11allb,subset=subpop))
  ggsc<-ggsurvplot_combine(ggls,xlab="OS Time (days)",title="Comparing the Primary Predictors Strength (MRI stage vs CRS): Cox (PH) Regression Model Comparisons by Survival Curves Estimate",legend.title = "Model by predictors", fin.df2[subpop,])#, palette = "Dark2")
  ggsc$plot<-ggsc$plot+scale_colour_manual(values = c(rgb(1,0,0,1),rgb(0,0,1,0.5),rgb(0.9,0.65,0)))
  return(ggsc$plot)
}
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

nostratmodel1<-coxph(ec.surv1~stage.MRI1+clinicrs+age1+tumour.vol1)
nostratmodel2<-coxph(ec.surv1~clinicrs+stage.MRI1+age1+tumour.vol1)
#anova(nostratmodel,nostratmodel1,test="LRT")
summary(stratmodel3)
summary(nostratmodel1)
summary(nostratmodel2)
summary(stratmodel3)
surv.df2<-data.frame(ec.surv1[,1],stage.MRI1,clinicrs,age1,tumour.vol1)
#library(Hmisc) #for rcorr
#library(sjPlot)
correlation::correlation(surv.df2,method="kendall")
nostratmodelwout1<-coxph(ec.surv1~clinicrs+age1+tumour.vol1)
nostratmodelwout2<-coxph(ec.surv1~stage.MRI1+age1+tumour.vol1)
anova(nostratmodel1,nostratmodelwout1,test="LRT")
anova(nostratmodel1,nostratmodelwout2,test="LRT")#signif.diff
AIC(nostratmodel1,nostratmodelwout2)
BIC(nostratmodel1,nostratmodelwout2)
summary(nostratmodel1)
summary(nostratmodelwout2)
#clinical rs is a better predictor among stage I patients

library(jtools)
library(ggstance)
library(broom.mixed)
cls=vector()
for(i in seq(5,70,by=5))
cls<-c(cls,colors(1)[i])
png("StageICoxModelComparison.png",width=16,height=8,units='in',res=300)
ps<-plot_summs(coxph.model1,coxph.model2,coxph.model3,coxph.model4,coxph.model5,model6,stratmodel,stratmodel2,stratmodel3,stratmodel4,model11,nostratmodel1,nostratmodelwout1,nostratmodelwout2,colors=cls,
           coefs = c("MRI Stage = IA"="stage.MRI1IA", "MRI Stage = IB"="stage.MRI1IB","MRI Stage = II"="stage.MRI1II","MRI Stage = IIIA"="stage.MRI1IIIA","MRI Stage = IIIB"="stage.MRI1IIIB","MRI Stage = IIIC1"="stage.MRI1IIIC1","MRI Stage = IIIC2"="stage.MRI1IIIC2","MRI Stage = IV"="stage.MRI1IV",
                     "Age"="age1", "Tumour Volume"="tumour.vol1",
                    "Tumour Grade = 1"="grade.histo11","Tumour Grade = 2"="grade.histo12","Tumour Grade = 3"="grade.histo13","Tumour Volume Grouping 1"="tumourv.groups1","Tumour Volume Grouping 2"="tumourv.groups2",
                    "Histological Type = endometrioid"="type.histo1endometrioid","Histological Type = clear cell"="type.histo1clear cell","Histological Type = mixed high grade"="type.histo1mixed high grade","Histological Type = serous"="type.histo1serous",
                    "Histological Stage = IA"="stage.histo1IA", "Histological Stage = IB"="stage.histo1IB","Histological Stage = II"="stage.histo1II","Histological Stage = IIIA"="stage.histo1IIIA","Histological Stage = IIIB"="stage.histo1IIIB","Histological Stage = IIIC1"="stage.histo1IIIC1","Histological Stage = IIIC2"="stage.histo1IIIC2","Histological Stage = IV"="stage.histo1IV",
                    "Clinical Risk Score = low"="clinicrslow","Clinical Risk Score = intermediate"="clinicrsintermediate","Clinical Risk Score = high"="clinicrshigh","Clinical Risk Score = advanced"="clinicrsadvanced",
                    "Tumour Volume > median (Grouping 2)"="tumourv.groups2>8839.888mm3",
                     "45-70 Age Group"="factor(age.groups1)45-70",">70 Age Group"="factor(age.groups1)>70"
           ))
ps+labs(title="Cox Models Comparison for Stage IA Subpopulation",subtitle ="based on the predictors common to at least two of the models")  #between the Predictors Common to at least Two Models")#overlapping predictors")
dev.off()

library(jtools)
library(ggstance)
library(broom.mixed)
#png("CoxModelComparisonW&WoutPrimaryPredictors.png",width=16,height=8,units='in',res=300)
model11b<-coxph(ec.surv1 ~ strata(stage.MRI1) + age1 + tumourv.groups2)
model11bc<-coxph(ec.surv1 ~ stage.MRI1 + age1 + tumourv.groups2)
png("CRSvsMRIstage-predictorstrength-CoxModelsComparisons.png",width=16,height=8,units='in',res=300)
ps<-plot_summs(model11nostrat,model11a,model11bc,colors=c("red","blue","orange"),coefs = c("MRI stage = IB relative to IA"="stage.MRI1IB",#"MRI Stage = II"="stage.MRI1II","MRI Stage = IIIA"="stage.MRI1IIIA","MRI Stage = IIIB"="stage.MRI1IIIB","MRI Stage = IIIC1"="stage.MRI1IIIC1","MRI Stage = IIIC2"="stage.MRI1IIIC2","MRI Stage = IV"="stage.MRI1IV",
                                                          "Age"="age1", "Tumour Volume"="tumour.vol1",
                                                          #"Tumour Grade = 1"="grade.histo11","Tumour Grade = 2"="grade.histo12","Tumour Grade = 3"="grade.histo13",
                                                          #"Tumour Volume Grouping 1"="tumourv.groups1",
                                                          #"Tumour Volume Grouping 2"="tumourv.groups2",
                                                          #"Histological Type = endometrioid"="type.histo1endometrioid","Histological Type = clear cell"="type.histo1clear cell","Histological Type = mixed high grade"="type.histo1mixed high grade","Histological Type = serous"="type.histo1serous",
                                                          #"Histological Stage = IA"="stage.histo1IA", "Histological Stage = IB"="stage.histo1IB","Histological Stage = II"="stage.histo1II","Histological Stage = IIIA"="stage.histo1IIIA","Histological Stage = IIIB"="stage.histo1IIIB","Histological Stage = IIIC1"="stage.histo1IIIC1","Histological Stage = IIIC2"="stage.histo1IIIC2","Histological Stage = IV"="stage.histo1IV",
                                                          "CRS = intermediate relative to low"="clinicrsintermediate","CRS = high relative to intermediate"="clinicrshigh","CRS = advanced relative to high"="clinicrsadvanced",
                                                          "Tumour Volume > (median = 8839.888 mm^3) relative to\ntumour volume <= 8839.888 mm^3 (grouping II)"="tumourv.groups2>8839.888mm3"
                                                          #"45-70 Age Group"="factor(age.groups1)45-70",">70 Age Group"="factor(age.groups1)>70"
),legend.title="Model by predictors",model.names = c("with MRI stage + CRS + age + tumour.groupingII","with CRS + age + tumour.groupingII","with MRI stage + age + tumour.groupingII"))
ps<-ps+labs(title="Comparing the Primary Predictors Strength (CRS vs MRI Stage):\nCox (PH) Regression Model Comparisons by Beta Coefficients",x="Estimate (Regression Coefficient)",y="Predictor")+theme(text=element_text(family="Arial"))#annotate(geom="text",label="n = 319, number of events = 57",x=Inf,y=-Inf)  #between the Predictors Common to at least Two Models")#overlapping predictors")
grid.arrange(ps,bottom="n = 319, number of events = 57 (subset: MRI stage = I)")
dev.off()

time1=ec.surv1[,1]
status=ec.surv1[,2]
datam11<-data.frame(time1,status,stage.MRI1,clinicrs,age1,tumourv.groups2)
summary(coxph(formula = Surv(time1,status) ~ stage.MRI1 + clinicrs + age1 + tumourv.groups2),
data=datam11)

stratmodel1<-coxph(ec.surv1~stage.MRI1+clinicrs+age1+tumourv.groups2)
stratmodelwout1<-coxph(ec.surv1~clinicrs+age1+tumourv.groups2)
stratmodelwout2<-coxph(ec.surv1~stage.MRI1+age1+tumourv.groups2)
plot_summs(stratmodel1,stratmodelwout1,stratmodelwout2,colors=cls,
           coefs = c("MRI Stage = IA"="stage.MRI1IA", "MRI Stage = IB"="stage.MRI1IB","MRI Stage = II"="stage.MRI1II","MRI Stage = IIIA"="stage.MRI1IIIA","MRI Stage = IIIB"="stage.MRI1IIIB","MRI Stage = IIIC1"="stage.MRI1IIIC1","MRI Stage = IIIC2"="stage.MRI1IIIC2","MRI Stage = IV"="stage.MRI1IV",
                     "Age"="age1", "Tumour Volume > median (Grouping 2)"="tumourv.groups2>8839.888mm3",
                     "Tumour Grade = 1"="grade.histo11","Tumour Grade = 2"="grade.histo12","Tumour Grade = 3"="grade.histo13",
                     #"Histological Type = endometrioid"="type.histo1endometrioid","Histological Type = clear cell"="type.histo1clear cell","Histological Type = mixed high grade"="type.histo1mixed high grade","Histological Type = serous"="type.histo1serous",
                     #"Histological Stage = IA"="stage.histo1IA", "Histological Stage = IB"="stage.histo1IB","Histological Stage = II"="stage.histo1II","Histological Stage = IIIA"="stage.histo1IIIA","Histological Stage = IIIB"="stage.histo1IIIB","Histological Stage = IIIC1"="stage.histo1IIIC1","Histological Stage = IIIC2"="stage.histo1IIIC2","Histological Stage = IV"="stage.histo1IV",
                     "Clinical Risk Score = low"="clinicrslow","Clinical Risk Score = intermediate"="clinicrsintermediate","Clinical Risk Score = high"="clinicrshigh","Clinical Risk Score = advanced"="clinicrsadvanced"
           ),legend.title="Model by primary predictors",model.names = c("with MRI Stage and CRS","without MRI Stage","without CRS"))
#CRS strongly predicts survival, stronger than MRI stage
summary(stratmodelwout1)#Concordance= 0.805  (se = 0.029 )
summary(stratmodelwout2)#Concordance= 0.73  (se = 0.034 )
summary(nostratmodelwout1)#Concordance= 0.81  (se = 0.029 )
summary(nostratmodelwout2)#Concordance= 0.744  (se = 0.034 )

model1<-coxph(ec.surv1~stage.MRI1+clinicrs)
model2<-coxph(ec.surv1~stage.MRI1+clinicrs+age1)
anova(model1,model2,test="LRT")#significantly diff
summary(model1)#Concordance= 0.74  (se = 0.036 )
summary(model2)#Concordance= 0.801  (se = 0.029 )
#this^suggests model2 is better
AIC(model1,model2)#lower for model2
BIC(model1,model2)#lower for model2
#they^suggest the same
model3<-coxph(ec.surv1~stage.MRI1+clinicrs+age1+tumour.vol1)
anova(model2,model3,test="LRT")#signif diff
summary(model2)#Concordance= 0.801  (se = 0.029 )
summary(model3)#Concordance= 0.813  (se = 0.029 )
AIC(model2,model3)#smaller for model3
BIC(model2,model3)#smaller for model3
#model3 is better
##
model4<-coxph(ec.surv1~stage.MRI1+clinicrs+age1+tumourv.groups2)
summary(model4)#xConcordance= 0.805  (se = 0.029 )
model10<-coxph(ec.surv1~stage.MRI1+clinicrs+strata(age1)+tumourv.groups2)#good!
summary(model10)
model11<-coxph(ec.surv1~strata(stage.MRI1)+clinicrs+age1+tumourv.groups2)
summary(model11)
summary(coxph(ec.surv1~stage.MRI1+strata(clinicrs)+age1+tumourv.groups2))

model5<-coxph(ec.surv1~stage.MRI1+clinicrs+age1+strata(tumourv.groups2))
summary(model5)#xConcordance= 0.791  (se = 0.03 )
model6<-coxph(ec.surv1~stage.MRI1+clinicrs+factor(age.groups1)+strata(tumourv.groups2))
summary(model6)#Preferred!Concordance= 0.769  (se = 0.031 )
##
model7<-coxph(ec.surv1~stage.MRI1+strata(clinicrs)+age1+tumour.vol1)
summary(model7)#X
model8<-coxph(ec.surv1~stage.MRI1+strata(clinicrs,age1)+tumour.vol1)
summary(model8)#X
model9<-coxph(ec.surv1~stage.MRI1+strata(clinicrs,age1)+tumourv.groups2)
summary(model9)#X!
summary(coxph(ec.surv1 ~ stage.MRI1 + strata(clinicrs)+ age1 + 
                tumourv.groups2))#X
summary(coxph(ec.surv1 ~ strata(stage.MRI1,clinicrs)+tumourv.groups2+age1))
#X


stageIBsubpop<-intersect(allInclStudy,which(corrected.data$X.2=="IB"))
length(stageIBsubpop)
ec.survb<-generate.Surv.obj(stageIBsubpop)
stage.MRIb<-corrected.data$X.2[stageIBsubpop]
ageb<-corrected.data$AgeAtDiagnosis[stageIBsubpop]
tumour.volb<-corrected.data$TumourVolume.mm.3.[stageIBsubpop]
tumourv.groupsb<-vector()
tumourv.groupsb[which(tumour.volb<=median(tumour.volb))]<-'<=8839.888mm3'
tumourv.groupsb[which(tumour.volb>median(tumour.volb))]<-'>8839.888mm3'
crsb<-corrected.data$Risk.score[stageIBsubpop]

stratmodelIB<-coxph(ec.survb~strata(stage.MRIb)+crsb+ageb+tumourv.groupsb)
summary(stratmodelIB)


ec.survall<-generate.Surv.obj(allInclStudy)
stage.MRIall<-corrected.data$X.2[allInclStudy]
table(stage.MRIall)
ageall<-corrected.data$AgeAtDiagnosis[allInclStudy]
tumour.volall<-corrected.data$TumourVolume.mm.3.[allInclStudy]
tumourv.groupsall<-vector()
tumourv.groupsall[which(tumour.volall<=median(tumour.volall))]<-'<=8839.888mm3'
tumourv.groupsall[which(tumour.volall>median(tumour.volall))]<-'>8839.888mm3'
crsall<-corrected.data$Risk.score[allInclStudy]

stratmodelall<-coxph(ec.survall~strata(stage.MRIall)+crsall+ageall+tumourv.groupsall)
summary(stratmodelall)
summary(coxph(ec.survall~stage.MRIall+crsall+ageall+tumour.volall))
model13<-coxph(ec.survall~stage.MRIall+crsall+strata(ageall)+tumour.volall)
summary(coxph(ec.survall~stage.MRIall+crsall+ageall+tumourv.groupsall))
model12<-coxph(ec.survall~stage.MRIall+crsall+strata(ageall)+tumourv.groupsall)
summary(model12)
summary(model13)
library(jtools)
library(ggstance)
library(broom.mixed)
png("StudyCoxModelComparison.png",width=16,height=8,units='in',res=300)
ps<-plot_summs(model12,model13,colors=cls,
               coefs = c("MRI Stage = IA"="stage.MRIallIA", "MRI Stage = IB"="stage.MRIallIB","MRI Stage = II"="stage.MRIallII","MRI Stage = IIIA"="stage.MRIallIIIA","MRI Stage = IIIB"="stage.MRIallIIIB","MRI Stage = IIIC1"="stage.MRIallIIIC1","MRI Stage = IIIC2"="stage.MRIallIIIC2","MRI Stage = IV"="stage.MRIallIV",
                         "Age"="ageall", "Tumour Volume"="tumour.volall",
                         "Tumour Grade = 1"="grade.histo11","Tumour Grade = 2"="grade.histo12","Tumour Grade = 3"="grade.histo13","Tumour Volume Grouping 1"="tumourv.groups1","Tumour Volume Grouping 2"="tumourv.groups2",
                         #"Histological Type = endometrioid"="type.histo1endometrioid","Histological Type = clear cell"="type.histo1clear cell","Histological Type = mixed high grade"="type.histo1mixed high grade","Histological Type = serous"="type.histo1serous",
                         #"Histological Stage = IA"="stage.histo1IA", "Histological Stage = IB"="stage.histo1IB","Histological Stage = II"="stage.histo1II","Histological Stage = IIIA"="stage.histo1IIIA","Histological Stage = IIIB"="stage.histo1IIIB","Histological Stage = IIIC1"="stage.histo1IIIC1","Histological Stage = IIIC2"="stage.histo1IIIC2","Histological Stage = IV"="stage.histo1IV",
                         "Clinical Risk Score = low"="crsalllow","Clinical Risk Score = intermediate"="crsallintermediate","Clinical Risk Score = high"="crsallhigh","Clinical Risk Score = advanced"="crsalladvanced",
                         "Tumour Volume > median (Grouping 2)"="tumourv.groups2>8839.888mm3"
                         
               ))
ps+labs(title="Comparison between Cox Models for the Entire Study Population",subtitle ="based on their common predictors")  #between the Predictors Common to at least Two Models")#overlapping predictors")
dev.off()


stagenIVsubpop<-setdiff(allInclStudy,which(corrected.data$X.2=="IV"))
length(stagenIVsubpop)
ec.survn4<-generate.Surv.obj(stagenIVsubpop)
stage.MRIn4<-corrected.data$X.2[stagenIVsubpop]
agen4<-corrected.data$AgeAtDiagnosis[stagenIVsubpop]
tumour.voln4<-corrected.data$TumourVolume.mm.3.[stagenIVsubpop]
tumourv.groupsn4<-vector()
tumourv.groupsn4[which(tumour.voln4<=median(tumour.voln4))]<-'<=8839.888mm3'
tumourv.groupsn4[which(tumour.voln4>median(tumour.voln4))]<-'>8839.888mm3'
crsn4<-corrected.data$Risk.score[stagenIVsubpop]

stratmodelnIV<-coxph(ec.survn4~strata(stage.MRIn4)+crsn4+agen4+tumourv.groupsn4)
summary(stratmodelnIV)

#coxph(ec.survb~stage.MRIb+ageb+strata(tumourv.groupsb)+clinicrsb)
# gga1<-ggadjustedcurves(nostratmodel1,variable = "stage.MRI1",data = surv.df2,reference = NULL,method = "average",fun = NULL)
# gga2<-ggadjustedcurves(nostratmodel1,variable = "clinicrs",data = surv.df2,reference = NULL,method = "average",fun = NULL)
# gga4<-ggadjustedcurves(nostratmodel2,variable = "clinicrs",data = surv.df2,reference = NULL,method = "average",fun = NULL)
# gga3<-ggadjustedcurves(nostratmodel2,variable = "stage.MRI1",data = surv.df2,reference = NULL,method = "average",fun = NULL)
# gga3a<-ggadjustedcurves(stratmodel3,variable = "stage.MRI1",data = surv.df2,reference = NULL,method = "average",fun = NULL)
# grid.arrange(gga1,gga3,gga3a,gga2,gga4, ncol = 3,nrow=2, top = "Compaison between primary explanatory variables in two models\n")

# png("ggplot_EC_OS_Cox-PHassumption412_updated.png",width=16,height=8,units='in',res=300)
# par(mfrow=c(3,1))
# #to see all plots at once (a plot for each predictor variable)
# plot(cox.zph(compl.model))
# dev.off()
# 
# png("ggplot_EC_OS_Cox-PHassumption412_whisto.png",width=16,height=8,units='in',res=300)
# par(mfrow=c(4,1))
# plot(cox.zph(coxph(Surv(os.timea,os.eventa)~stage.MRI+agep+tumour.vol+grade.histo)))
# dev.off()

# png("ggplot_EC_OS_Cox-PHassumption412_normalized.png",width=16,height=8,units='in',res=300)
# par(mfrow=c(3,1))
# plot(cox.zph(coxph(Surv(os.timea,os.eventa)~stage.MRI.rescalef+age.rescale+tumour.vol.rescale)))
# dev.off()

#normalized<-cox.zph(coxph(Surv(os.timea,os.eventa)~stage.MRI.rescalef+age.rescale+tumour.vol.rescale))
#png("ggplot_EC_OS_Cox-PHassumption412_normalized_line_thru0.png",width=16,height=8,units='in',res=300)
#par(mfrow=c(3,1))

# par(mfrow=c(1,1))
# plot(normalized[1])
# abline(h=0,col="red")
# 
# par(mfrow=c(1,1))
# plot(normalized[2])
# abline(h=0,col="red")
# 
# par(mfrow=c(1,1))
# plot(normalized[3])
# abline(h=0,col="red")


# par(mfrow=c(1,1))
# plot(cox.zph(compl.model)[1])
# abline(h=0,col="red")
# 
# par(mfrow=c(1,1))
# plot(cox.zph(compl.model)[2])
# abline(h=0,col="red")
# 
# par(mfrow=c(1,1))
# plot(cox.zph(compl.model)[3])
# abline(h=0,col="red")

#abline(h=0,col=2)#col from color
#abline(h=0,col=3)
#the dash lines are confidence intervals/confidence bands around that smoother line
#add a (red) line at 0 to see how often is a change of 0 contained in this interval 
#dev.off()

fit.Cox<-function(){
  os.timef<-os.time1[-213]
  os.eventf<-os.event1[-213]
  library(survival)
  ec.osf<-Surv(os.timef,os.eventf)
  MRI.stagef<-data$X.2[allincl]
  km_stage <- survfit(ec.osf ~ MRI.stagef)
  plot(km_stage, fun = "cloglog", xlab = "Time (in days) using log",
       ylab = "log-log survival", main = "log-log curves by clinic") 
  library(survminer)
  ggsurvplot(km_stage, fun = "cloglog", data=data)
  
  histo.stagef<-data$STAGE[allincl]
  km_histo.stage<-survfit(ec.osf ~ histo.stagef)
  plot(km_histo.stage, fun = "cloglog", xlab = "Time (in days) using log",
       ylab = "log-log survival", main = "log-log curves by clinic") 
  #convergent and divergent lines -> important covariate(s) missed
  #boxplot(ec.osf~)
  plot(survfit(ec.osf ~ histo.gradef), fun = "cloglog", xlab = "Time (in days) using log",
       ylab = "log-log survival", main = "log-log curves by clinic") 
  ggsurvplot(survfit(ec.osf ~ histo.gradef), fun = "cloglog", data=data)
  ggsurvplot(survfit(ec.osf ~ histo.typef), fun = "cloglog", data=data)
  
  
  library(survival)
  #coxph fits the Cox (PH) regression model 
  coxph.stage<-coxph(ec.osf ~ MRI.stagef,data=data)
  summary(coxph.stage)
  exp(confint(coxph.stage))
  test.ph_stage <- cox.zph(coxph.stage)
  ind3b<-intersect(which(data$X.2=="IIIB"),allincl)
  ind3c<-intersect(which(data$X.2=="IIIC"),allincl)
  #no events in the subset of cases with a positive value of stage IIIB and IIIC, respectively
  #beta values too small (high in module) for them
  data$Death.Date[ind3b]
  data$Death.Date[ind3c]
  
  agef<-data$AgeAtDiagnosis[allincl]
  
  histo.gradef<-data$Grade[allincl]
  histo.typef<-data$X.1[allincl]
  coxphmodel<-coxph(ec.osf ~ agef+histo.gradef+histo.typef)
  test.ph <- cox.zph(coxphmodel)
  #plot(test.ph)
  ggcoxzph(test.ph)
  coxph_histo.gradef<-coxph(ec.osf~histo.gradef)
  coxph_histo.typef<-coxph(ec.osf~histo.typef)
  cox.zph(coxph_histo.gradef)
  cox.zph(coxph_histo.typef)
  
  #here we can assume the proportional hazards, since p-vals are not
  #statistically significant
  library(survminer)
  ggcoxzph(test.ph)#here we can observe that beta does not vary much over time
  ggcoxdiagnostics(coxphmodel, type = "dfbeta",
                   linear.predictions = FALSE, ggtheme = theme_bw())
  ggcoxdiagnostics(coxphmodel, type = "deviance",
                   linear.predictions = FALSE, ggtheme = theme_bw())
  ggcoxfunctional(ec.osf ~ agef + log(agef) + sqrt(agef),data=data)
  #it only works for categorical covariates
  # km_age <- survfit(ec.osf ~ agef)
  # plot(km_age, fun = "cloglog", xlab = "Time (in days) using log",
  #      ylab = "log-log survival", main = "log-log curves by clinic") 
  # library(survminer)
  # ggsurvplot(km_age, fun = "cloglog", data=data)
  deaths<-length(which(os.eventf==1))
  #78 deaths & 3-5 covariates -> 15-26 events per covariate (so at least 10)
}

add.missing.deathInfo<-function()
{
  colnames(revisedDeathInfo)[1]<-"PatientID"
  updatedPatients<-revisedDeathInfo$PatientID[which(revisedDeathInfo$X=="not dead")]
  corrected.data$PatientID[corrected.data$PatientID %in% updatedPatients]
  corrected.data$Death.Date[corrected.data$PatientID %in% updatedPatients]
  corrected.data$Death.Date[corrected.data$PatientID %in% updatedPatients]<-"no"
  corrected.data$Death.Date[which(corrected.data$Death.Date=="No")]<-"no"
  table(corrected.data$Death.Date)
  #corrected.data[20,2],corrected.data[20,53]
  ids<-which(corrected.data$PatientID %in% updatedPatients)
  library(xlsx)
  workbook_vmar <- loadWorkbook(file = "EC_data_10mar2021.xlsx")
  sheets <- getSheets(workbook_vmar)
  ids3<-ids+3
  rows  <- getRows(sheets$AllCases,rowIndex = ids3)   # get all the rows
  
  cc<-getCells(rows,colIndex = c(53))
    #setCellValue(cc$`627.53`,"hey")
    #nm<-names(cc)[70]
    #setCellValue(cc[[nm]],"hey")
  
  indc <- paste0(ids+3,".53")
  # col1<-paste("\`627",".53\`",sep="")
  # cc$as.name(col1)
  # cc1<-cc
  # lapply(cc[[nms]],setCellValue,"no")
  
  
  new.indc<-which(names(cc) %in% indc)
  #names(cc)<-indc
  for(i in new.indc){
    nmi<-names(cc)[i]
    setCellValue(cc[[nmi]],"no")
  }
  #setCellValue(cc[[60280]],"no")
  # idp<-4
  # cc1[paste0(ids,".53")]<-rep(c(no),times=length(ids))
  # setCellValue(cc[paste0(ids+3,".53")],rep(c(no),times=length(ids)))
  # mapply(setCellValue, cc[paste0(ids,".53")], no)
  # lapply(cc[paste0(ids,".53")],setCellValue,no)
  # for(cell in cc[paste0(ids,".53")])
  #   print(cell)
  cc[paste0(ids+3,".53")]<-rep(c(no),times=length(ids))
  # cells<-cc[paste0(ids+3,".53")]
  # setCellValue(cc,cc)
  # lapply(cells, setCellValue, "no")
  # indc <- paste0(ids+3,".53")
  # mapply(setCellValue, cc[indc],"no")
  # lapply(cc[indc],function(x){})
  # for(i in indc)
  #   setCellValue(cc[i],"no")
  # lapply(setCellValue,lapply(cc,function(x){setCellValue(x,no)}))
    #setCellValue(cell,no)
  #cc1[c("14.30","14.32")]<-rep(c(no),2)
  # lapply(cc1,function(x){
  #   vname<-paste(idp,".53",sep="")
  #   cc[[vname]]<-ifelse(idp %in% ids,no,x)
  #   idp<-idp+1
  # })
  #tm<-paste0(627,".53")
  #options(useFancyQuotes = c("\x60","\x60","\x60","\x60"))#"\xab", "\xbb", "\xbf", "?"))
  #col1<-as.name(sQuote(tm))#,dQuote(".53"))
  
  no<-new(J("java.lang.String"), "nope")
  cc[setdiff(indc,names(cc))]<-rep(c(no),times=length(ids)-3)
  setCellValue(cc[["135.53"]],"no")
  cc$"135.53"<-"no"
  saveWorkbook(workbook_vmar,file = "EC_data_12mar2021.xlsx")
  #how many empty cells are in Death Date column
  #emptydd<-which(data$Reasons.for.Exclusion %in% c("include","missingMRI")&data$Death.Date=="")
  #length(emptydd)
  studypop1<-studypop
  studypop<-which(corrected.data$Reasons.for.Exclusion=="include" &
          (grepl("0",corrected.data$Death.Date) | corrected.data$Death.Date %in% c("no","No")))#,"no ")))
  #not in version 15 march
  corrected.data$Risk.score[65]
  corrected.data$Reasons.for.Exclusion[65]<-"no risk score assessment"
  #corrected.data$X.2[2]<-"IIIC1"
  table(corrected.data$Reasons.for.Exclusion[which(corrected.data$X.1 %in% excluded.subtypes)])
  corrected.data$PatientID[intersect(which(corrected.data$Reasons.for.Exclusion=="missingMRI"),which(corrected.data$X.1 %in% excluded.subtypes))]
  corrected.data$Reasons.for.Exclusion[178]<-"missingMRI & stromal"
  #not in version 15 march
  
  table(corrected.data$Grade[studypop])
  table(corrected.data$STAGE[studypop])
  table(corrected.data$X.1[studypop])
  table(corrected.data$X.2[studypop])
  table(corrected.data$X.2)
  allstudypop1<-allstudypop
  allstudypop<-which(corrected.data$X.2!="" & corrected.data$Reasons.for.Exclusion %in% c("include","missingMRI") &
                    (grepl("0",corrected.data$Death.Date) | corrected.data$Death.Date %in% c("no","No")))#,"no ")))
  table(corrected.data$Grade[allstudypop])
  table(corrected.data$STAGE[allstudypop])
  table(corrected.data$X.1[allstudypop])
  table(corrected.data$X.2[allstudypop])
  table(corrected.data$Availability.of.MRI[allstudypop])
  table(corrected.data$Reasons.for.Exclusion[allstudypop])
  table(chemotherapy[allstudypop])#in column "Yes.No"
  table(radiotherapy[allstudypop])#in column "Yes.No.1"
  
  genstage<-unique(c(2,intersect(allstudypop,which(corrected.data$X.2=="IIIC")),intersect(allstudypop,which(corrected.data$STAGE=="IIIC"))))
  allcids<-c(grep("PatientID",colnames(corrected.data)),grep("radiomicsCaseID",colnames(corrected.data)),grep("AgeAtDiagnosis",colnames(corrected.data)),grep("OpDate",colnames(corrected.data)),grep("STAGE",colnames(corrected.data)),grep("Grade",colnames(corrected.data)),grep("X.1",colnames(corrected.data)),grep("MRI.date",colnames(corrected.data)),cids)
  genst.df<-data.frame()
  genst.df<-corrected.data[genstage,allcids]
  create.ExcelTable.file(genst.df,"general_stageIIIC_insteadOf_IIIC1_or_IIIC2")
  corrected.data$Reasons.for.Exclusion[genstage]
  corrected.data$Availability.of.MRI[genstage]
  #not in version 15 March
  2 %in% genstage
  corrected.data$X.2[2]<-"IIIC"
  corrected.data$Reasons.for.Exclusion[2]<-"incorrect staging classification"
  corrected.data$STAGE[26]<-"IV"
  corrected.data$X.2[26]<-"IIIC1"
  corrected.data$Reasons.for.Exclusion[26]<-"include"
  corrected.data$STAGE[36]
  corrected.data$X.2[36]<-"IIIC1"
  corrected.data$Reasons.for.Exclusion[36]<-"include"
  corrected.data$STAGE[170]<-"IIIC2"
  corrected.data$X.2[170]
  corrected.data$Reasons.for.Exclusion[170]
  corrected.data$STAGE[208]
  corrected.data$X.2[208]
  corrected.data$Reasons.for.Exclusion[208]
  
  
  chemotherapy<-corrected.data$"Yes.No"
  y.n<-corrected.data$"Yes.No"
  chemotherapy[which(grepl("y",y.n,ignore.case = TRUE)&!grepl("\\?",y.n)&nchar(y.n)<=nchar("y (for recurrence)"))]<-"yes"
  chemotherapy[which(grepl("n",y.n,ignore.case = TRUE)&!grepl("\\?",y.n)&nchar(y.n)<=3|grepl("decline",y.n))]<-"no"
  table(chemotherapy)
  chemo.complete<-which(chemotherapy %in% c("yes","no"))
  missingChemo<-setdiff(studypop,chemo.complete)
  radiotherapy<-corrected.data$"Yes.No.1"
  y.n.1<-corrected.data$"Yes.No.1"
  radiotherapy[which(grepl("y",y.n.1,ignore.case = TRUE)&!grepl("\\?",y.n.1)&nchar(y.n.1)<=nchar("y(Mt Vernon)"))]<-"yes"
  radiotherapy[which((grepl("n",y.n.1,ignore.case = TRUE)&!grepl("\\?",y.n.1)&nchar(y.n.1)<=3|grepl("no |n ",y.n.1))&!grepl("y/n",radiotherapy))]<-"no"
  table(radiotherapy)
  radth.complete<-which(radiotherapy %in% c("yes","no"))
  missingRadth<-setdiff(studypop,radth.complete)
  length(studypop)-length(missingRadth)
  table(corrected.data$ImageVoxelSize.mm.)
}

find.citations<-function(){
  citation("survival")
  citation("survminer")
}

add.missing.MRInfo<-function(){
  pid1<-data71MRIstaged$PatientID[which(is.na(data71MRIstaged$MRI.Stage))]
  ids1<-which(corrected.data$PatientID %in% pid1)
  corrected.data$X.2[ids1]<-NA
  corrected.data$Availability.of.MRI[ids1]
  corrected.data$Reasons.for.Exclusion[ids1]
  corrected.data$Availability.of.MRI[469]<-"missing MRI stage ONLY"
  corrected.data$Reasons.for.Exclusion[469]<-"missingMRI"
  #manually changed (forgot about it at first)
  corrected.data$Availability.of.MRI[179]<-"missing MRI stage, depth of invasion, and other MRI findings"
  
  library(xlsx)
  workbook_vmar1<-loadWorkbook(file = "EC_data_11mar2021.xlsx")
  sheets1<-getSheets(workbook_vmar1)
  Ids1<-ids1+3
  rows1<-getRows(sheets1$AllCases,rowIndex = Ids1)
  cc1<-getCells(rows1,colIndex = c(126))
  for(i in 1:length(cc1)){
    nmi<-names(cc1)[i]
    setCellValue(cc1[[nmi]],"NA",showNA = TRUE)
  }
  saveWorkbook(workbook_vmar1,file = "EC_data_11mar2021_ARexcl.xlsx")
  
  inclMRIs<-which(data71MRIstaged$Myometrial.Invasion %in% c("sup","deep"))
  pid2<-data71MRIstaged$PatientID[inclMRIs]
  ids2<-which(corrected.data$PatientID %in% pid2)
  #corrected.data$Depth.of.Myometrial.Invasion[ids2]<-data71MRIstaged$Myometrial.Invasion[inclMRIs]
  s<-grep("MRI.Stage", colnames(data71MRIstaged))
  e<-grep("L.inguinal.LN", colnames(data71MRIstaged))
  start<-grep("Y.N.4", colnames(corrected.data))
  end<-grep("Y.N.22", colnames(corrected.data))
  data71MRIstaged$MRI.Stage[which(data71MRIstaged$MRI.Stage=="1A")]<-"IA"
  for(i in s:e)
    print(table(data71MRIstaged[inclMRIs,i]))
  cids<-c(grep("X.2",colnames(corrected.data)),grep("Depth.of.Myometrial.Invasion",colnames(corrected.data)),seq(start,end,by=5))
  corrected.data[ids2,cids]<-data71MRIstaged[inclMRIs,seq(s,e)]
  corrected.data$Availability.of.MRI[ids2]<-"yes"
  corrected.data$Reasons.for.Exclusion[ids2]<-"include"
  cids2<-c(cids,grep("Availability.of.MRI",colnames(corrected.data)),grep("Reasons.for.Exclusion",colnames(corrected.data)))
  
  library(xlsx)
  workbook_vmar2<-loadWorkbook(file = "EC_data_11mar2021_ARexcl.xlsx")
  sheets2<-getSheets(workbook_vmar2)
  Ids2<-ids2+3
  rows2<-getRows(sheets2$AllCases,rowIndex = Ids2)
  cc2<-getCells(rows2,colIndex = cids2)
  library(stringr)
  for(i in 1:length(cc2)){
    nmi<-names(cc2)[i]
    vs<-str_extract_all(nmi, "[[:digit:]]+")
    rnb<-as.numeric(vs[[1]][1])-3
    cnb<-as.numeric(vs[[1]][2])
    setCellValue(cc2[[nmi]],corrected.data[rnb,cnb])
  }
  saveWorkbook(workbook_vmar2,file = "EC_data_11mar2021_ARincl.xlsx")

  nomyom<-which(grepl("none",data71MRIstaged$Myometrial.Invasion))
  pid3<-data71MRIstaged$PatientID[nomyom]
  ids3<-which(corrected.data$PatientID %in% pid3)
  corrected.data$Availability.of.MRI[ids3]
  corrected.data$Reasons.for.Exclusion[ids3]
  corrected.data$Depth.of.Myometrial.Invasion[ids3]
  data71MRIstaged$Myometrial.Invasion[which(data71MRIstaged$Myometrial.Invasion=="none ")]<-"none"
  for(i in s:e)
    print(table(data71MRIstaged[nomyom,i]))
  for(i in cids)
    print(table(corrected.data[ids3,i]))
  corrected.data[ids3,cids]<-data71MRIstaged[nomyom,seq(s,e)]
  
  library(xlsx)
  workbook_vmar3<-loadWorkbook(file = "EC_data_11mar2021_ARincl.xlsx")
  sheets3<-getSheets(workbook_vmar3)
  Ids3<-ids3+3
  rows3<-getRows(sheets3$AllCases,rowIndex = Ids3)
  cc3<-getCells(rows3,colIndex = cids)
  library(stringr)
  for(i in 1:length(cc3)){
    nmi<-names(cc3)[i]
    vs<-str_extract_all(nmi, "[[:digit:]]+")
    rnb<-as.numeric(vs[[1]][1])-3
    cnb<-as.numeric(vs[[1]][2])
    setCellValue(cc3[[nmi]],corrected.data[rnb,cnb])
  }
  saveWorkbook(workbook_vmar3,file = "EC_data_11mar2021_ARfin.xlsx")

  corrected.data$Reasons.for.Exclusion[510]<-"missingMRI & breast cancer"
  #pe positiile ids3 si ultima din ids1 avem none
  #renal lesions
  #2 identical patients
  
  
  exclMRIs<-which(grepl("EXCLUDE",data40MRIstage$MRI.Stage))
  data40MRIstaged$MRI.Stage[exclMRIs]
  data40MRIstaged$Reasons.for.Exclusion[exclMRIs]
  table(corrected.data$Depth.of.Myometrial.Invasion)
  data40MRIstaged$Radiomics.CaseID[seq(37,76)]
  data40MRIstaged$Myometrial.Invasion[which(data40MRIstaged$Myometrial.Invasion=="<50%")]<-"sup"
  data40MRIstaged$Myometrial.Invasion[which(data40MRIstaged$Myometrial.Invasion==">50%")]<-"deep"
  toBfilled<-setdiff(seq(37,76),exclMRIs)
  pid4<-data40MRIstaged$Radiomics.CaseID[toBfilled]
  ids4<-which(corrected.data$radiomicsCaseID %in% pid4)
  cids4<-c(grep("MRI.date",colnames(corrected.data)),cids)
  fill.ids<-c(grep("MRI.Date",colnames(data40MRIstaged)),seq(s,e))
  corrected.filled<-toBfilled[order(match(data40MRIstaged$Radiomics.CaseID[toBfilled], corrected.data$radiomicsCaseID[ids4]))]
  corrected.data[ids4,cids4]<-data40MRIstaged[corrected.filled,fill.ids]
  
  library(xlsx)
  workbook_vmar4<-loadWorkbook(file = "EC_data_11mar2021_ARfin.xlsx")
  sheets4<-getSheets(workbook_vmar4)
  Ids4<-ids4+3
  rows4<-getRows(sheets4$AllCases,rowIndex = Ids4)
  cc4<-getCells(rows4,colIndex = cids4)
  library(stringr)
  for(i in 1:length(cc4)){
    nmi<-names(cc4)[i]
    vs<-str_extract_all(nmi, "[[:digit:]]+")
    rnb<-as.numeric(vs[[1]][1])-3
    cnb<-as.numeric(vs[[1]][2])
    setCellValue(cc4[[nmi]],corrected.data[rnb,cnb])
  }
  saveWorkbook(workbook_vmar4,file = "EC_data_11mar2021_NBlatest.xlsx")
}

extract.out.revisedMRIs<-function(){
  rindices<-seq(37,76)
  radiomics.ids<-data40MRIstaged$Radiomics.CaseID[rindices]
  tbextracted<-which(dataMissedMRI$Radiomics.CaseID %in% radiomics.ids)
  new.dataMissedMRI<-dataMissedMRI[-tbextracted,]
  library(writexl)
  write_xlsx(new.dataMissedMRI,"C:\\Users\\Casi\\Documents\\EndoProj\\hadMRI_but_missingMRIdata_updated.xlsx")
}

recheck.MRI.avail<-function()
{
  table(data$MRI.date)
  table(data$Custom.worklist)
  which(data$Availability.of.MRI=="yes")
  
  #discrepancies between the scan type and MRI date
  
  #there is an MRI date, but it states CT only
  MRIdate_noMRIscan<-which(data$MRI.date!=""&!grepl("MRI",data$Custom.worklist)&grepl("CT",data$Custom.worklist))
  #there is no MRI date, but it states that MRI has been performed
  MRIscan_noMRIdate<-which(data$MRI.date==""&grepl("MRI",data$Custom.worklist))
  ii<-c(MRIdate_noMRIscan,MRIscan_noMRIdate)
  
  PatientID<-data$PatientID[ii]
  Included.Excluded<-data$Reasons.for.Exclusion[ii]
  Radiomics.CaseID<-data$radiomicsCaseID[ii]
  AgeAtDiagnosis<-data$AgeAtDiagnosis[ii]
  MRI.Date<-data$MRI.date[ii]
  Custom.Worklist<-data$Custom.worklist[ii]
  MRI.Availability<-data$Availability.of.MRI[ii]
  MRI.Stage<-data$X.2[ii]
  MRI.myom.inv<-MRI.minv[ii]
  #ds=date/scan discrepancies
  df.ds<-data.frame(PatientID,Included.Excluded,Radiomics.CaseID,AgeAtDiagnosis,MRI.Date,Custom.Worklist,MRI.Availability,MRI.Stage)
  df.ds$"Myometrial Invasion"<-MRI.myom.inv
  df.ds
  start<-grep("Y.N.4", colnames(data))
  end<-grep("Y.N.22", colnames(data))
  cn<-c("Other involvement: Cervical stroma","Vagina/Parametria","Adnexa","Serosa","Peritoneum/Omentum","R-pelvic LN","L-pelvic LN","Para-aortic LN","R-inguinal LN","L-inguinal LN")
  length(cn)
  df.ds[cn]<-data[ii,seq(start,end,by=5)]
  df.ds
  # "MRI date is recorded, but only CT scan is mentioned"
  # "missing MRi date and data, but it is mentioned that patient had MRI"
  # "missing MRI date, but evidence of MRI is recorded"
  df.ds$"Error"<-
  c(rep("no MRI scan mentioned even though MRI date is given",length(MRIdate_noMRIscan)),
    rep("no MRI date but MRI scan is mentioned",length(MRIscan_noMRIdate)))

  library(writexl)
  write_xlsx(df.ds,"C:\\Users\\Casi\\Documents\\EndoProj\\MRI_date_vs_scanType.xlsx")
}

find.missing.stage.with.MRI<-function(){
  MRInostage_findings<-which(data$Availability.of.MRI=="missing MRI stage ONLY")
  #there is MRI, but MRI stage is missing
  MRInostage_date<-which(data$MRI.date!="" & data$X.2=="") 
  MRInostage_any_info<-which(grepl("missing",data$Availability.of.MRI))
  
  #patient had MRI but the MRI staging info is missing
  MRInostage<-which(grepl("missing",data$Availability.of.MRI)&data$MRI.date!="")
  data$Availability.of.MRI[MRInostage]
  df.dst<-create.sep.df(MRInostage)
  df.dst
  df.dst$"Error"<-c("MRI date recorded but MRI staging information is missing")

  library(writexl)
  write_xlsx(df.dst,"C:\\Users\\Casi\\Documents\\EndoProj\\MRIdate_noMRIstaging.xlsx")  
}

which(data$Availability.of.MRI=="missing y/n in MRI findings")
which(dataMissedMRI$Reasons.for.Exclusion!="missingMRI")

create.sep.df<-function(ii){
  PatientID<-data$PatientID[ii]
  Included.Excluded<-data$Reasons.for.Exclusion[ii]
  Radiomics.CaseID<-data$radiomicsCaseID[ii]
  AgeAtDiagnosis<-data$AgeAtDiagnosis[ii]
  MRI.Date<-data$MRI.date[ii]
  Custom.Worklist<-data$Custom.worklist[ii]
  MRI.Availability<-data$Availability.of.MRI[ii]
  MRI.Stage<-data$X.2[ii]
  MRI.myom.inv<-MRI.minv[ii]
  #ds=date/scan discrepancies
  df<-data.frame(PatientID,Included.Excluded,Radiomics.CaseID,AgeAtDiagnosis,MRI.Date,Custom.Worklist,MRI.Availability,MRI.Stage)
  df$"Myometrial Invasion"<-MRI.myom.inv
  
  start<-grep("Y.N.4", colnames(data))
  end<-grep("Y.N.22", colnames(data))
  cn<-c("Other involvement: Cervical stroma","Vagina/Parametria","Adnexa","Serosa","Peritoneum/Omentum","R-pelvic LN","L-pelvic LN","Para-aortic LN","R-inguinal LN","L-inguinal LN")
  df[cn]<-data[ii,seq(start,end,by=5)]
  return(df)
}

add.column.ExcelTable<-function(col.name){
  library(xlsx)
  wb<-loadWorkbook("EC_data_14feb2021_provera.xlsx")
  shh<-getSheets(wb)
  col.to.find<-paste0("^",col.name,"$")
  col.indx<-grep(col.to.find, colnames(data))
  addDataFrame(data[,col.indx,drop=F], shh$AllCases, startColumn = col.indx, startRow = 3, row.names=FALSE)
  saveWorkbook(wb,"EC_data_15feb2021.xlsx")
}


#l<-list()
l<-vector(mode = "list", length = nrow(data))
filterCa.out<-function(cname){
  # otherCa<-apply(data, 2, 
  #                function(x) 
  #                { 
  #                  x[str_extract_all(x, cname)!="character(0)"]
  #                })
  # vec<-vector()
  # for(i in 1:length(otherCa))
  # {
  #   if(length(otherCa[[i]])==1 && !is.na(otherCa[[i]]))
  #     vec<-c(vec,(otherCa[[i]]))
  # }
  # return(vec)
  
  vec<-vector()
  for(i in 1:nrow(data))
  {
    cell.val<-data[i,][grepl(cname,data[i,],ignore.case=TRUE)]
    # we could have used [data[,i] %like% "renal"] for case sensitive
    if(length(cell.val)>=1){
      vec<-c(vec,cell.val)
      l[[i]]<<-cell.val # global 
    }
    else
      l[[i]]<-NULL
  }
  return(vec)
}

l3<-vector(mode = "list", length = nrow(data))
genetic.syndromes<-function()
{
  lynch<-c(filterCa.out("Lynch"),filterCa.out("HNPCC"))
  Filter(Negate(is.null), l)
  
  l2<-l
  l2[which(lapply(l,is.null) == F)]<-"Lynch syndrome"
  #l2
  
  l<<-vector(mode = "list", length = nrow(data))
  
  brca2<-filterCa.out("BRCA")
  Filter(Negate(is.null), l)
  
  l3<<-l2
  names(l3)<-names(l2)
  l3[which(lapply(l2,is.null) == T&lapply(l,is.null) == F)]<<-"BRCA2 gene mutation"
  l3[[188]]
  #l3
  
  b<-l2[which(lapply(l2,is.null) == F&lapply(l,is.null) == F)]
  if(length(b)>0){
    b<-paste(b,"& BRCA2 gene mutation", sep=" ")
    # add the new cancer in the latest list
    l3[which(lapply(l2,is.null) == F&lapply(l,is.null) == F)]<<-b
    #l3[[188]]
  }
  #b
}

#add.categ.genetics
#insert.end.col<-function(colname,inFile,outFile,lname)
#{
  data$"Genetic Predispositions"<- lapply(l3,function(x){ifelse(is.null(x),"none",x)})
  data$"Genetic Predispositions"#colname, lname
  library(xlsx)
  workbook <- loadWorkbook("updated_EC_data - Copy.xlsx")#"EC_cleaned_anonymised_data_latest - Copy.xlsx"
  sheets <- getSheets(workbook)
  sheets
  #columnToPreserve = readColumns(sheets$AllCases, 16, 16, startRow=1)
  #addDataFrame(columnToPreserve, sheets$test, startColumn = 17, row.names = FALSE)
  #and so on...
  #ncol(data)
  #colnames(data)
  #simulating a dataframe with a single column
  #ncol(data)+1
  addDataFrame(data[,131,drop=F], sheets$AllCases, startColumn = 131, startRow = 3, row.names=FALSE)
  #write.xlsx(df, workbook, sheetName="AllCases")#, row.names=TRUE)
  saveWorkbook(workbook, 'EC_data_28jan2021.xlsx')
  
#}

init.vars()
clean.stage()
clean.minv()
write.csv(accuracy.data,"stagingVSinvasion.csv",row.names=TRUE)
print.incorrect.stage()
print.incorrect.minv()
stageIdata<-subset(accuracy.data,histo.stage %in% c("IA","IB") & MRI.stage %in% c("IA","IB"))
stageIdata
number.stage.stI<-sum(stageIdata$histo.stage!=stageIdata$MRI.stage)
number.stage.stI
number.minv.stI<-sum(stageIdata$histo.minv!=stageIdata$MRI.minv)
number.minv.stI
index.stage.stI<-which(stageIdata$histo.stage!=stageIdata$MRI.stage)
index.minv.stI<-which(stageIdata$histo.minv!=stageIdata$MRI.minv)
index.stage.stI
index.minv.stI
#add.column(stageIdata,"Yes.No.1","radiotherapy")


renalCa<-filterCa.out("renal")
renalCa
renal.cancer<-renalCa[!grepl("adrenal",renalCa) & !grepl("infrarenal",renalCa) & !grepl("renal lesion",renalCa)]
# adrenal, infrarenal, renal lesion should be excluded
renal.cancer

ind<-0
l<<-lapply(l,function(x){ ind<<-ind+1
  if(ind<=length(l)){
    if(!str_contains(l[[ind]],"adrenal")&!str_contains(l[[ind]],"infrarenal")&!str_contains(l[[ind]],"renal lesion"))
    { 
      x<<-l[[ind]] 
      #lfin[[ind]]<-"renal cancer"
    }
    else
    {
      x<<-NULL
      #lfin[[ind]]<-NULL
    }
  }
})
names(l)<<-names(renal.cancer)
Filter(Negate(is.null), l)
l

lfin<-vector(mode = "list", length = length(l))
#lfin<-l
#lfin<-list()
lfin[which(lapply(l,is.null) == F)]<-"renal cancer"
lfin
# ind<<-0
# lfin<-lapply(lfin,function(x){ ind<<-ind+1
#   if(!is.null(x))
#     lfin[[ind]]<-"renal cancer"
#   else
#     lfin[[ind]]<-NULL
# })
#names(lfin)<<-names(renal.cancer)
Filter(Negate(is.null), lfin)
lfin

ovarianCa<-filterCa.out("ovar")
ovarianCa
ovary.cancer<-ovarianCa[grepl("ovarian",ovarianCa,ignore.case=TRUE)&!grepl("ovarian vein",ovarianCa)&!grepl("ovarian mets",ovarianCa) | grepl("tumour",ovarianCa,ignore.case=TRUE) | ovarianCa=="ovary ca"]
#ovarian vein; ovary ca (particular cases)
ovary.cancer

ind<<-0
l<<-lapply(l,function(x){ ind<<-ind+1
if(ind<=length(l)){
  if(str_contains(l[[ind]],"ovarian",ignore.case = TRUE)&!str_contains(l[[ind]],"ovarian vein")&!str_contains(l[[ind]],"ovarian mets") | str_contains(l[[ind]],"tumour", ignore.case = TRUE) | str_contains(l[[ind]],"ovary ca"))
  {
      x<<-l[[ind]] 
      print(x)
  }
  else
  {
      x<<-NULL
  }
}
})
names(l)<<-names(ovary.cancer)
Filter(Negate(is.null), l)

lfin
l
lfin2<-lfin
lfin2
lfin2[which(lapply(lfin,is.null) == T&lapply(l,is.null) == F)]<-"ovarian cancer"
lfin2

length(lfin2)
b<-lfin[which(lapply(lfin,is.null) == F&lapply(l,is.null) == F)]
if(length(b)>0){
  b<-paste(b,"& ovarian cancer", sep=" ")
  lfin2[which(lapply(lfin,is.null) == F&lapply(l,is.null) == F)]<-b
}
b


breastCa<-filterCa.out("breast")
breastCa # no need to filter
names(l)<<-names(breastCa)
l
lfin3<-lfin2
lfin3[which(lapply(lfin2,is.null) == T&lapply(l,is.null) == F)]<-"breast cancer"
lfin3

b<-lfin2[which(lapply(lfin2,is.null) == F&lapply(l,is.null) == F)]
if(length(b)>0){
  b<-paste(b,"& breast cancer", sep=" ")
  lfin3[which(lapply(lfin2,is.null) == F&lapply(l,is.null) == F)]<-b
}
b


lungCa<-filterCa.out("lung")
lungCa
#lung(s) met(s) /lesion, "lobectomy lung", mets lung, "pleura, lungs" 
#"lymph nodes and lung", "lung and peritoneum", "lung and liver", "lung"
#"lung, peritoneal spleen", "lung + pelvic side wall", "umbilical nodule and lung"
#lung nodule(s), progression lungs
#recurrence
lung.cancer<-lungCa[grepl("lung primary",lungCa,ignore.case=TRUE) | grepl("lung ca",lungCa)]
lung.cancer

ind<<-0
l<<-lapply(l,function(x){ ind<<-ind+1
if(ind<=length(l)){
  if(str_contains(l[[ind]],"lung primary",ignore.case = TRUE) | str_contains(l[[ind]],"lung ca"))
  {
    x<<-l[[ind]] 
    print(x)
  }
  else
  {
    x<<-NULL
  }
}
})
#names(l)<<-names(lung.cancer)
Filter(Negate(is.null), l)
#they seem like 2 different entries, but they have the same index

lfin4<-lfin3
lfin4[which(lapply(lfin3,is.null) == T&lapply(l,is.null) == F)]<-"lung cancer"
lfin4

b<-lfin3[which(lapply(lfin3,is.null) == F&lapply(l,is.null) == F)]
if(length(b)>0){
  b<-paste(b,"& lung cancer", sep=" ")
  lfin4[which(lapply(lfin3,is.null) == F&lapply(l,is.null) == F)]<-b
}
b


colonCa<-filterCa.out("colon cancer")
colon.cancer<-c(colonCa,filterCa.out("concurrent sigmoid primary"))
colon.cancer
# l<-vector(mode = "list", length = nrow(data))
# for(i in 1:nrow(data))
# {
#   cell.val1<-data[i,][grepl("colon cancer",data[i,])]
#   cell.val2<-data[i,][grepl("concurrent sigmoid primary",data[i,])]
#   if(length(cell.val1)>=1){
#     l[[i]]<-cell.val1 
#   }
#   else
#     l[[i]]<-NULL
#   if(length(cell.val2)>=1){
#     l[[i]]<-cell.val2
#     print(cell.val2)
#   }
#   else
#     l[[i]]<-NULL
# }
#names(l)<<-names(colon.cancer)
Filter(Negate(is.null), l)

lfin5<-lfin4
lfin5[which(lapply(lfin4,is.null) == T&lapply(l,is.null) == F)]<-"colon cancer"
lfin5[[411]]

b<-lfin4[which(lapply(lfin4,is.null) == F&lapply(l,is.null) == F)]
if(length(b)>0){
  b<-paste(b,"& colon cancer", sep=" ")
  lfin5[which(lapply(lfin4,is.null) == F&lapply(l,is.null) == F)]<-b
}
b


colorectalCa<-filterCa.out("colorectal")
Filter(Negate(is.null), l)

lfin6<-lfin5
lfin6[which(lapply(lfin5,is.null) == T&lapply(l,is.null) == F)]<-"colorectal cancer"
lfin6[[181]]

b<-lfin5[which(lapply(lfin5,is.null) == F&lapply(l,is.null) == F)]
if(length(b)>0){
  b<-paste(b,"& colorectal cancer", sep=" ")
  lfin6[which(lapply(lfin5,is.null) == F&lapply(l,is.null) == F)]<-b
}
b


analCa<-filterCa.out("anal")
Filter(Negate(is.null), l)

lfin7<-lfin6
lfin7[which(lapply(lfin6,is.null) == T&lapply(l,is.null) == F)]<-"anal cancer"
lfin7[[181]]
l[[181]]

b<-lfin6[which(lapply(lfin6,is.null) == F&lapply(l,is.null) == F)]
if(length(b)>0){
  b<-paste(b,"& anal cancer", sep=" ")
  lfin7[which(lapply(lfin6,is.null) == F&lapply(l,is.null) == F)]<-b
  lfin7
}
b


lynch<-c(filterCa.out("Lynch"),filterCa.out("HNPCC"))
Filter(Negate(is.null), l)

lfin8<-lfin7
lfin8[which(lapply(lfin7,is.null) == T&lapply(l,is.null) == F)]<-"Lynch syndrome"
lfin8

# see if our current filter points out an entry 
#that is already marked with some cancers;
# extract which cancers have been detected previously
b<-lfin7[which(lapply(lfin7,is.null) == F&lapply(l,is.null) == F)]
if(length(b)>0){
  b<-paste(b,"& Lynch syndrome", sep=" ")
  # add the new cancer in the latest list
  lfin8[which(lapply(lfin7,is.null) == F&lapply(l,is.null) == F)]<-b
  lfin8[[51]]
  lfin8[[181]]
}
b


brca2<-filterCa.out("BRCA")
Filter(Negate(is.null), l)

lfin9<-lfin8
lfin9[which(lapply(lfin8,is.null) == T&lapply(l,is.null) == F)]<-"BRCA2 gene mutation"
lfin9[[188]]
lfin9

b<-lfin8[which(lapply(lfin8,is.null) == F&lapply(l,is.null) == F)]
if(length(b)>0){
  b<-paste(b,"& BRCA2 gene mutation", sep=" ")
  # add the new cancer in the latest list
  lfin9[which(lapply(lfin8,is.null) == F&lapply(l,is.null) == F)]<-b
  lfin9[[188]]
}
b


grade<-data$"X.1"
grade
lfin10<-lfin9
length(lfin10)
excluded.subtypes<-c("no cancer","uterine metastasis from other cancer","stromal","sarcoma","no histology","typical HP","atypical HP","neuroendocrine")
for(nECsubtype in excluded.subtypes){
  for(idx in which(grade==nECsubtype))
    if(is.null(lfin10[[idx]]))
      lfin10[[idx]]<-nECsubtype
    else{
      b<-lfin10[[idx]]
      b<-paste(b,"&",nECsubtype, sep=" ")
      lfin10[[idx]]<-b
      print(idx)
    }
}
which(grade=="uterine metastasis from other cancer")
length(lfin10)
lfin10[[218]]
lfin10[[438]]
lfin10[[443]]
lfin10[[392]]
lfin10


data$"Reasons for exclusion"<- lapply(lfin10,function(x){ifelse(is.null(x),"include",x)})
data$"Reasons for exclusion"
library(xlsx)
workbook <- loadWorkbook("EC_cleaned_anonymised_data_latest - Copy.xlsx")
sheets <- getSheets(workbook)
sheets
#columnToPreserve = readColumns(sheets$AllCases, 16, 16, startRow=1)
#addDataFrame(columnToPreserve, sheets$test, startColumn = 17, row.names = FALSE)
#and so on...
ncol(data)
colnames(data)
#simulating a dataframe with a single column
addDataFrame(data[,130,drop=F], sheets$AllCases, startColumn = 130, startRow = 3, row.names=FALSE)
#write.xlsx(df, workbook, sheetName="AllCases")#, row.names=TRUE)
#saveWorkbook(workbook, 'updated_EC_data.xlsx')
#saveWorkbook(workbook, 'updated_EC_data - Copy.xlsx')

genetic.syndromes()
l3[[442]]#507,441,442 - Lynch syndrome
l3[[188]]#BRCA gene mutation 


#class(colnames(data)[1])
#data$Scan.CT[data$Scan.CT %like% "yes - recurrence"]

# shopping_list <- c("apples x4", "bag of flour", "bag of sugar", "milk x2")
# shopping_list[str_extract_all(shopping_list,"r")!="character(0)"]
# class(shopping_list)
# shopping_list[length(str_extract_all(shopping_list,"r"))==1]
# str_extract_all(shopping_list,"r")[1]=="character(0)"
#select(data,"X.3")

