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

pie_preproc<-function(ss_pie_raw=NULL,filter_freechoice=T,only_firstfree=F,usemeanprior=F){
  print(unique(ss_pie_raw$ID))
  numseg<-max(ss_pie_raw$num_segments)
  ss_pie_scon<-split(ss_pie_raw,ss_pie_raw$con_num)
  indexsx<-rbind(data.frame(segnum=1:numseg,type="samplehx"),
                 #data.frame(segnum=1:numseg,type="rewardhx"),
                 data.frame(segnum=1:numseg,type="v_bayes"),
                 data.frame(segnum=1:numseg,type="choice"))
  indexsx$variname<-paste0(indexsx$type,indexsx$segnum)
  commenvir<-as.environment(list())
  assign("vbayvars",indexsx$variname[indexsx$type=="v_bayes"],envir = commenvir)
  assign("samplehxvars",indexsx$variname[indexsx$type=="samplehx"],envir = commenvir)
  assign("choicevars",indexsx$variname[indexsx$type=="choice"],envir = commenvir)
  
  ss_proc<-do.call(rbind,lapply(ss_pie_scon,function(sx){
    #Okay REDESIGN!!!!
    sxw<-sx[which(sx$RT!=0),]
    
    ext<-lapply(sxw$trial,function(i){
      #print(i)
      storaget<-as.environment(list())
      
      segchoice<-sxw[i,"selected_segment"]
      segrwad<-sxw[i,"win"]
      
      #Current Choice;
      choicearray<-rep(0,numseg)
      choicearray[segchoice]<-1
      assign("choicearray",choicearray,envir = storaget)
      #Past Sample History;
      assign("samplehxarray",sapply(1:numseg,function(xj){
        if(i!=1){
          return(length(which(sxw[1:(i-1),"selected_segment"]==xj)))}else{
            return(0)
          }
      }),envir = storaget)
      #Value calculated by bayes rule
      assign("vbayarray",sapply(1:numseg, function(y) {
        nreward<-sum(sxw[1:i,"win"])
        nchoice<-length(which(sxw[1:i,"selected_segment"]==y))
        nchoicegivenrewar<-length(which(sxw[1:i,"selected_segment"]==y & sxw[1:i,"win"]==1))
        ntotal<-i
        pcgivenreward<-ifelse(nreward==0,0,(nchoicegivenrewar / nreward))
        if(usemeanprior){
          preward<-mean(sxw$selected_prob)} else {preward<-nreward/ntotal}
        pchoice<-(nchoice/ntotal)
        if(pchoice!=0){
          v_bayes <- (pcgivenreward * preward) / pchoice
        }else{v_bayes <- 0}
        return(v_bayes)}),envir = storaget)
      
      ext_df<-do.call(cbind,lapply(c("vbay","samplehx","choice"),function(jx) {
        arrayx<-as.data.frame(as.list(get(paste0(jx,"array"),envir = storaget)),col.names = get(paste0(jx,"vars"),envir = commenvir))
        selectedx<-as.data.frame(as.list(get(paste0(jx,"array"),envir = storaget)[segchoice]),col.names = paste0(jx,"_selected"))
        return(cbind(arrayx,selectedx))
      }))
      return(ext_df)
    }) 
    sxw<-cbind(sxw,do.call(rbind,ext))
    if(max(sxw$num_segments) < numseg) {sxw[indexsx$variname[indexsx$segnum > max(sxw$num_segments)]]<-NA}
    if(only_firstfree) {
      sxw$firstfree<-FALSE
      sxw$firstfree[max(sxw$num_segments)+1]<-TRUE
    }
    
    sxw$u_l<-1-(sxw$samplehx_selected/(sxw$trial-1))
    # for (ix in c(indexsx$variname,"samphx","rewhx")) {
    #   d<-sxw[,ix]
    #   sxw[paste0(ix,"_lag")]<-dplyr::lag(d)
    #   sxw[paste0(ix,"_lead")]<-dplyr::lead(d)
    # }
    #print(length(sxw))
    return(sxw)
  }))
  
  if(filter_freechoice){ss_proc<-ss_proc[!as.logical(ss_proc$forced_choice),]}
  if(only_firstfree) {ss_proc<-ss_proc[ss_proc$fristfree,]}
  return(ss_proc)
}
