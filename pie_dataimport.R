#Behavioral Processing for Pie Task:
source("pie_utility.R")


if(file.exists("pie_data.rdata")) {load("pie_data.rdata")} else {
boxsyncpath<-"/Volumes/bek/Box Sync"
boxsyncpath<-findbox()
piedata_raw<-pie_getdata(boxsyncpath)
}

pie_data_proc<-ProcApply(piedata_raw$list,pie_preproc,filter_freechoice=F,only_firstfree=F)

pie_data_proc_f<-ProcApply(piedata_raw$list,pie_preproc,filter_freechoice=T,only_firstfree=F)

pie_firstfree<-ProcApply(piedata_raw$list,pie_preproc,filter_freechoice=T,only_firstfree=T)$df

save(piedata_raw,pie_data_proc,pie_data_proc_f,pie_firstfree,file = "pie_data.rdata")
