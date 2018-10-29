#Pie Utility Functions:

#General Functions
getifswitched<-function(y) {
  c(y[1],y[1:length(y)-1])->y_lag
  return(!y==y_lag)
}
#Label probability function
lableVar<-function(dfx) {
  if (length(grep("if",names(dfx)))>0) {
    for (jx in grep("if*",names(dfx))) {
      as.logical(dfx[[jx]])->temp
      dfx[which(temp),jx]<-gsub("if","",names(dfx)[jx])
      dfx[which(!temp),jx]<-paste0("Not_",gsub("if","",names(dfx)[jx]))
    } }
  return(dfx)
}
#clean up list function
cleanuplist<-function(listx){
  if (any(sapply(listx, is.null))){
    listx[sapply(listx, is.null)] <- NULL}
  return(listx)
}
#Generate probability function
genProbability<-function(dfx,condition=c("Context","Emotion"),response=c("FaceResponseText"),excludeNA=T,missresp=NA) {
  if (excludeNA) {
    if (is.na(missresp)) {
      dfx<-dfx[which(!is.na(dfx[[response]])),] } else {dfx<-dfx[which(dfx[[response]]!=missresp),]}
  }
  dfx<-droplevels(dfx)
  #whichone<-c(condition,response)
  interaction(dfx[condition])->interactions
  
  nwx<-do.call(rbind,lapply(attributes(as.factor(dfx[[response]]))$levels, function(resp) {
    prob<-data.frame(
      p=sapply(attributes(interactions)$levels, function(x) {
        ( length(which(as.character(dfx[[response]])==resp & interactions==x)) / length(which(interactions==x)) ) -> px
        return(px)
      }),
      resp=resp)
    for (n in 1:length(condition)) {
      prob[condition[n]]<-sapply(strsplit(rownames(prob),split = ".",fixed = T),"[[",n)
    }
    rownames(prob)<-NULL
    prob$ID<-unique(dfx$ID)
    lableVar(prob)
  }) )
  
  return(nwx)
}
#sigma transformation
sigmatransform<-function(x) {
  (1./(1+exp(x)))->y
  return(y)
}
#Get the depth of a list
depthoflist <- function(list,thisdepth=0){
  if(!is.list(list)){
    return(thisdepth)
  }else{
    return(max(unlist(lapply(list,depthoflist,thisdepth=thisdepth+1))))    
  }
}
#Add centerscale to a list
addcenterscaletolist<-function(list) {
  test<-lapply(list, scale,center=T,scale=T)
  names(test)<-paste(names(list),"centerscaled",sep = "_")
  newlist<-c(list,test) 
  return(newlist)
}
#VIF function for lme
vif.lme <- function (fit) {
  ## adapted from rms::vif
  v <- vcov(fit)
  nam <- names(fixef(fit))
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)] }
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v }
#Find Box
findbox<-function() {
  if (Sys.getenv("USER")=="jiazhouchen") {boxdir <- "/Users/jiazhouchen/Box Sync"
  } else if (Sys.getenv("USER")=="jiazhou") {boxdir <- "/Volumes/bek/Box Sync"} else {
    boxdir<-system("find ~ -iname 'Box*' -maxdepth 2 -type d",intern = T)}
  return(boxdir)
}

ProcApply<-function(listx=NULL,FUNC=NULL,...) {
  proc_af<-lapply(X = listx,FUN = FUNC,... = ...)
  return(list(list=proc_af,
         df=do.call(rbind,proc_af)))
}


#Pie Specific;
pie_getdata<-function(boxsyncpath=NULL){
  pieroot<-file.path(boxsyncpath,"skinner","data","matlab task data","pie_task")
  piedata_raw<-lapply(list.files(path = pieroot,pattern = ".*_outstruct.csv",full.names = T),read.csv)
  names(piedata_raw)<-gsub("([0-9]+).*$", "\\1",list.files(path = pieroot,pattern = ".*_outstruct.csv"))
  
  piedata_raw<-lapply(gsub("([0-9]+).*$", "\\1",list.files(path = pieroot,pattern = ".*_outstruct.csv")),function(ID){
    rawdata<-read.csv(file.path(pieroot,paste0(ID,"_outstruct.csv")))
    rawdata$ID<-ID
    return(rawdata)
  })
  
  piedata_raw_all<-do.call(rbind,piedata_raw)
  
  return(list(list=piedata_raw,df=piedata_raw_all))
}

pie_preproc<-function(ss_pie_raw=NULL,filter_freechoice=T,only_firstfree=F){
  numseg<-max(ss_pie_raw$num_segments)
  ss_pie_scon<-split(ss_pie_raw,ss_pie_raw$con_num)
  ss_proc<-do.call(rbind,lapply(ss_pie_scon,function(sx){
    indexsx<-rbind(data.frame(segnum=1:numseg,type="samplehx"),
                   data.frame(segnum=1:numseg,type="rewardhx"),
                   data.frame(segnum=1:numseg,type="choice"))
    indexsx$variname<-paste0(indexsx$type,indexsx$segnum)
    tw<-as.data.frame(as.list(rep(0,3*numseg)))
    names(tw)<-indexsx$variname
    sxw<-merge(sx,tw,all = T)
    for (i in sx$trial) {
      segchoice<-sxw[i,"selected_segment"]
      segrwad<-sxw[i,"win"]
      samplevar<-indexsx$variname[indexsx$type=="samplehx" & indexsx$segnum==segchoice]
      rewvar<-indexsx$variname[indexsx$type=="rewardhx" & indexsx$segnum==segchoice]
      choicevar<-indexsx$variname[indexsx$type=="choice" & indexsx$segnum==segchoice]
      if (i==1) {
        choice_hx<-0
        rew_hx<-0
      } else {
        choice_hx<-sxw[(i-1),samplevar]
        rew_hx<-sxw[(i-1),rewvar]
      }
      sampleupdate<-choice_hx+1
      rew_update<-(rew_hx+segrwad)/sampleupdate
      sxw[(i),samplevar]<-sampleupdate
      sxw[(i),rewvar]<-rew_update
      sxw[i:length(sxw[[1]]),indexsx$variname]<-sxw[i,indexsx$variname]
      sxw[(i),choicevar]<-1
    }
    if(max(sxw$num_segments) < numseg) {sxw[indexsx$variname[indexsx$segnum > max(sxw$num_segments)]]<-NA}
    if(only_firstfree) {
      sxw$fristfree<-FALSE
      sxw$fristfree[max(sxw$num_segments)+1]<-TRUE
    }
    for (ix in indexsx$variname) {
      d<-sxw[,ix]
      sxw[paste0(ix,"_lag")]<-dplyr::lag(d)
      sxw[paste0(ix,"_lead")]<-dplyr::lead(d)
    }
    return(sxw)
  }))
  
  if(filter_freechoice){ss_proc<-ss_proc[!as.logical(ss_proc$forced_choice),]}
  if(only_firstfree) {ss_proc<-ss_proc[ss_proc$fristfree,]}
  return(ss_proc)
}
