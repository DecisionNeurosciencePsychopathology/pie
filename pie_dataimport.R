#Behavioral Processing for Pie Task:
source("pie_utility.R")
#Step 1, data importation:
boxsyncpath<-"/Volumes/bek/Box Sync"

piedata_raw<-pie_getdata(boxsyncpath)

pie_data_proc<-ProcApply(piedata_raw$list,pie_preproc,filter_freechoice=T,only_firstfree=F)

pie_firstfree<-ProcApply(piedata_raw$list,pie_preproc,filter_freechoice=T,only_firstfree=T)$df




