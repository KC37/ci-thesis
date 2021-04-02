#header2 <- scan("EC_anonymised_data_latest - Copy.csv", nlines = 2, what = character())
# header 3
#data<-read.csv("EC_anonymised_data_latest - Copy.csv",header=TRUE,skip=2)
data<-read.csv("EC_data_22feb2021.csv",header=TRUE,skip=2)
head(data)
colnames(data)

corrected.data<-read.csv("EC_data_10mar2021.csv",header=TRUE,skip=2)
colnames(corrected.data)

dataS2<-read.csv("EC_data_15feb2021Sheet2.csv",header=TRUE)
colnames(dataS2)

dataMissedMRI<-read.csv("hadMRI_but_missingMRIdata.csv",skip=1,header=TRUE)
colnames(dataMissedMRI)

dataEC416<-read.csv("EC416WseparatedCRcateg&RRcateg.csv",header=TRUE)
colnames(dataEC416)

data30MRI<-read.csv("First 30 Rockall_AR_for casiana.csv",header=TRUE)
colnames(data30MRI)

data40MRIstaged<-read.csv("hadMRI_but_missingMRIstage_NB_FOR CASIANA.csv",skip=1,header=TRUE)
colnames(data40MRIstaged)

data71MRIstaged<-read.csv("hadMRI_but_missingMRIdata_updated_AR.csv",skip=1,header=TRUE)
colnames(data71MRIstaged)

revisedDeathInfo<-read.csv("D_edit_Extra_Patients_RadiomicsStudy-144from416-forDiana.csv",header=TRUE)
colnames(revisedDeathInfo)

source("stageIdepth.R")
