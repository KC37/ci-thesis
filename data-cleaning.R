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

# comparison between stages IA (superficial) and IB (deep) and other stages against their corresponding involvement
#incorrectly staging based on the MRI findings according to FIGO table 2009
compareMRIstageW.MRIfindings<-function(){
  # compare depth of invasion against stage I: IA vs IB
  
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

         
compareHistoW.risk.scores<-function(){
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

         
check4discrep<-function()
{
  compareMRIstageW.MRIfindings()
  compareHistoW.risk.scores()
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
  }

         
compare.myTable.wExtTable<-function(){
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
  #a column from an external table with time from operation to 3 august 2020(=the censorship date) even for the patients 
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
}
  
         
handle.patientsheet2<-function(){  
  #incomplete!! added a complete query later on
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

search4genetic.pred<-function(){
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
  
}
         
exhaustic.search4other.malignancies<-function(){
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
}

