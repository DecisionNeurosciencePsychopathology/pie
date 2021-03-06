#Behavioral Processing for Pie Task:
#setwd("~/code/pie")

source("pie_utility.R")
overwrite_pie_data=1
dataset='PSU'


if(file.exists("pie_data.rdata")&&overwrite_pie_data==0) {load("pie_data.rdata")} else {
  #boxsyncpath<-"C:/Users/brownvm2/Box Sync" #/Volumes/bek/Box Sync"
  #boxsyncpath<-"~/Box Sync" #/Volumes/bek/Box Sync"
  #boxsyncpath<-findbox()
  boxsyncpath='C:/Users/brownvm2/Box'
  piedata_raw<-pie_getdata(dataset,boxsyncpath)
}
pie_data_proc<-ProcApply(multicorenum = 1,piedata_raw$list,pie_preproc,filter_freechoice=F,only_firstfree=F,usemeanprior=F)

#pie_data_proc_f<-ProcApply(piedata_raw$list,pie_preproc,filter_freechoice=T,only_firstfree=F)

pie_firstfree<-ProcApply(multicorenum = 1,piedata_raw$list,pie_preproc,filter_freechoice=T,only_firstfree=T)$df

varyingvars<-names(pie_data_proc$df)[grepl("[1-8]",names(pie_data_proc$df)) & !names(pie_data_proc$df) %in% names(piedata_raw$df)]
pie_data_proc_long<-reshape(data = pie_data_proc$df,v.names = unique(gsub("[0-9.]", "", varyingvars)),
        varying = varyingvars, idvar = c("ID","trial","block_num"),times = 1:8,timevar = "segment",
        direction = "long")
pie_data_proc_long<-pie_data_proc_long[order(pie_data_proc_long$ID,pie_data_proc_long$block_num,pie_data_proc_long$trial),]
pie_data_proc_long<-pie_data_proc_long[as.logical(apply(pie_data_proc_long[unique(gsub("[0-9.]", "", varyingvars))],1,function(x) { any(!is.na(x)) })),]

save(piedata_raw,pie_data_proc,pie_firstfree,pie_data_proc_long,file = "pie_data.rdata")
