#Behavioral Processing for Pie Task:
setwd("~/code/pie")

source("pie_utility.R")


if(file.exists("pie_data.rdata")) {load("pie_data.rdata")} else {
boxsyncpath<-"/Volumes/bek/Box Sync"
boxsyncpath<-findbox()
piedata_raw<-pie_getdata(boxsyncpath)
}
pie_data_proc<-ProcApply(multicorenum = 8,piedata_raw$list,pie_preproc,filter_freechoice=F,only_firstfree=F,usemeanprior=F)

#pie_data_proc_f<-ProcApply(piedata_raw$list,pie_preproc,filter_freechoice=T,only_firstfree=F)

pie_firstfree<-ProcApply(piedata_raw$list,pie_preproc,filter_freechoice=T,only_firstfree=T)$df

save(piedata_raw,pie_data_proc,pie_firstfree,file = "pie_data.rdata")
