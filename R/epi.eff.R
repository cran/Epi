## Program to calculate effects
## Michael Hills August 2006
## version 0.1.3

epi.eff<-function(response,type="metric",exposure,strata=NULL,
                  control=NULL,fup=NULL,data=NULL) {


  ## usage

  if(missing(response)&missing(exposure)&missing(strata)&
      missing(control)&missing(fup)&missing(data))  {
      cat("Usage:","\n")
      cat("\n")
      cat("epi.eff (","\n")
      cat("response =              ,","\n")
      cat("type     = \"metric\"     ,","   (one of metric/binary/failure/count)","\n")
      cat("exposure =              ,","\n")
      cat("strata   =              ,","   (must be a factor)","\n")
      cat("control  = data.frame() ,","\n")
      cat("fup      =              ,","   (must be numeric)","\n")
      cat("data     =              ,","   (must be a data frame)","\n")
      cat(")","\n")
      cat("\n")
      stop("No arguments specified")
  }
  
  ## attaches the dataframe specified in data=
  
  if(!is.null(data)) {
    attach(data,2)
    on.exit(detach(pos=2))
  }

  ## stores the variable names for response, etc.

  rname <-deparse(substitute(response))
  tname<-deparse(substitute(type))
  ename<-deparse(substitute(exposure))
  sname<-deparse(substitute(strata))

  ## performs a few checks

  if(!is.numeric(response))stop("Response must be numeric, not a factor")
  if(rname==ename)stop("Same variable specified as response and exposure")
  if(rname==sname)stop("Same variable specified as response and strata")
  if(sname==ename)stop("Same variable specified as strata and exposure")
  
  if(missing(response))stop("Must specify the response","\n")
  if(missing(exposure))stop("Must specify the exposure","\n")
  if(!is.null(strata)&!is.factor(strata))stop("Stratifying
    variable must be a factor")
  if(type=="binary") {
    tmp<-(response==0 | response==1)
    if(all(tmp,na.rm=TRUE)==FALSE)
      stop("Binary response must be coded 0,1 or NA")
  }
  if(type=="failure") {
    tmp<-(response==0 | response==1)
    if(all(tmp,na.rm=TRUE)==FALSE)
      stop("Failure response must be coded 0,1 or NA")
  }
  if(type=="failure"&missing(fup))stop("Must specify a follow-up variable
    when type is failure")
    
  
  
  ## prints out some information about variables

  cat("***************************************************************","\n")
  cat("response      : ", rname, "\n")
  cat("type          : ", tname, "\n")
  cat("exposure      : ", ename, "\n")
  if(!is.null(control))cat("control vars  : ",names(control),"\n")
  if(!is.null(strata)) {
    cat("stratified by : ",sname,"\n")
  }
  cat("\n")
  if(is.factor(exposure)) {
    cat(ename,"is a factor with levels: ")
    cat(levels(exposure),"\n")
  }
  else {
    cat(ename,"is metric","\n")
  }
  if(!is.null(strata)) {
    cat(sname,"is a factor with levels: ")
    cat(levels(strata),"\n")
  }
  cat("\n")
  if(type=="metric")cat("effects are measured as differences in means","\n")
  if(type=="binary")cat("effects are measured as odds ratios","\n")
  if(type=="failure")cat("effects are measured as rate ratios","\n")
# cat("***************************************************************","\n")
# cat("\n")

  ## translates type of response into family

  if (type=="metric") family<-"gaussian"
  if (type=="binary") family<-"binomial"
  if (type=="failure") family<-"poisson"

  ## gets number of levels for exposure if a factor

  if(is.factor(exposure)) {
    nlevE<-length(levels(exposure))
  }

  ## labels the output
  
  if(is.factor(exposure)) {
    cat("effect of",ename,"on",rname,"\n")
  }
  else {
    cat("effect of an increase of 1 unit in",ename,"on",rname,"\n")
  }
  if(!is.null(control)) {
    cat("controlled for",names(control),"\n")
  }
  if(!is.null(strata)) {
    cat("stratified by",sname,"\n")
  }

  ## metric or binary response
  
  if (type=="metric" | type=="binary") {

  ## no stratifying variable
    
      if(is.null(strata)) {
          if(is.null(control)) {            
            m<-glm(response~exposure,family=family)
            mm<-glm(response~1,family=family,subset=!is.na(exposure))
          }
          else  {
              m<-glm(response~.+exposure,family=family,
                     subset=!is.na(exposure),data=control)
              mm<-glm(response~.,family=family,
                      subset=!is.na(exposure),data=control)
          }
          if (type=="metric") {
              res<-ci.lin(m,subset=c("Intercept","exposure"))
              res<-res[,c(1,5,6)]
            }
          if(type=="binary") {
              res<-ci.lin(m,subset=c("Intercept","exposure"),Exp=TRUE)
              res<-res[,c(5,6,7)]
            }
          res<-signif(res,3)
          colnames(res)<-c("Effect","2.5%","97.5%")
          if(is.factor(exposure)) {
            rownames(res)[2:nlevE]<-paste("level",2:nlevE,"vs 1")
          }

          aov <- anova(mm,m,test="Chisq")
          return(list(res,paste("Test for no effects of exposure on",
                 aov[2,3],"df:","p=",format.pval(aov[2,5],digits=3))))
      }
      
   ## stratifying variable
      
      if(!is.null(strata)) {
        nlevS<-length(levels(strata))
          if(is.null(control)) {
            m<-glm(response~strata/exposure,family=family)
            mm<-glm(response~strata+exposure,family=family)
          }
          else {
            m <-glm(response~strata/exposure + .,family=family,
            data=control)
            mm <-glm(response~strata+exposure + .,family=family,
            data=control)
          }
          if(type=="binary") {
            res<-ci.lin(m,subset=c("strata"),Exp=TRUE)[c(-1:-(nlevS-1)),
                        c(5,6,7)]
          }
          else {
              res<-ci.lin(m,subset=c("strata"))[c(-1:-(nlevS-1)),c(1,5,6)]
          }
          res<-signif(res,3)
          colnames(res)<-c("Effect","2.5%","97.5%")
          if(is.factor(exposure)) {
            newrownames<-NULL
            for(i in c(1:(nlevE-1))) {
              newrownames<-c(newrownames,
                             paste("strata",1:nlevS,"(level",i+1,"vs 1)"))
            }
          }
          else {
            newrownames<-paste("strata",1:nlevS)
          }
          rownames(res)<-newrownames
          aov<-anova(mm,m,test="Chisq")
          return(list(res,paste("Test for effect modification on",
          aov[2,3],"df:","p=",format.pval(aov[2,5],digits=3))))
       }
  }

## failure or count response
  
  if (type=="failure" | type=="count") {

    ## no strata
    
      if (is.null(strata)) {
          if (is.null(control)) {
              m<-glm(response~exposure+offset(log(y)),family=family)
              mm<-glm(response~1+offset(log(y)),family=family,
                      subset=!is.na(exposure))
          }
          else  {
              m<-glm(response~.+exposure+offset(log(y)),family=family,
              data=control)
              mm<-glm(response~.+offset(log(y)),family=family,
                      subset=!is.na(exposure),data=control)
          }
              res<-ci.lin(m,subset=c("Intercept","exposure")
                          ,Exp=TRUE)[,c(5,6,7)]
              res<-signif(res,3)
              colnames(res)<-c("Effect","2.5%","97.5%")
              if(is.factor(exposure)) {
                nlevE<-length(levels(exposure))
                rownames(res)[2:nlevE]<-paste("level",2:nlevE,"vs 1")
              }
              aov<-anova(mm,m,test="Chisq")
          return(list(res,paste("Test for no effects of exposure on",
          aov[2,3],"df:","p=",format.pval(aov[2,5],digits=3))))
      }

      ## strata specified
      
      if(!is.null(strata)) {
        nlevS<-length(levels(strata))
          if(is.null(control)) {
            m<-glm(response~strata/exposure+offset(log(y)),
                   family=family)
            mm<-glm(response~strata+exposure+offset(log(y)),
                    family=family)
          }
          else {
            m <-glm(response~.+strata/exposure+offset(log(y)),
                    family=family,data=control)
            mm<-glm(response~.+strata+exposure+offset(log(y)),
                    family=family,data=control)
          }
          res<-ci.lin(m,subset=c("strata"),Exp=TRUE
          )[c(-1:-(nlevS-1)),c(5,6,7)]
          res<-signif(res,3)
          colnames(res)<-c("Effect","2.5%","97.5%")
          if(is.factor(exposure)) {
            newrownames<-NULL
            for(i in c(1:(nlevE-1))) {
              newrownames<-c(newrownames,
                             paste("strata",1:nlevS,"(level",i+1,"vs 1)"))
            }
          }
          else {
            newrownames<-paste("strata",1:nlevS)
          }
          rownames(res)<-newrownames
          aov<-anova(mm,m,test="Chisq")
          return(list(res,paste("Test for effect modification on",
          aov[2,3],"df:","p=",format.pval(aov[2,5],digits=3))))
       }
  }
}

