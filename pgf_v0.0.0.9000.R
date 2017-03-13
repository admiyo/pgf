#################
# load support files
source(file="ipak.R")
source(file="SL_Wrapper.R")
e1<-read.table("Rpgf_dat.txt", sep="\t",header=T)


# load packages
packages <- c("SuperLearner","data.table","doParallel","doRNG","sampling","plyr","profvis")
ipak(packages)

# create SL library
SL.library<-c("SL.mean","SL.glm","SL.nnet",
              "SL.gam.1","SL.gam.2","SL.gam.3","SL.gam.4",
              "SL.glmnet","SL.glmnet.0.25","SL.glmnet.0.5","SL.glmnet.0.75")
#################

# set monte carlo re-sample size
montecarlo<-1228

# arrange dataset & create variables
data_format<-function(data, identifier, timescale){
  m<-data
  m$id<-m[,identifier]
  m$mm<-m[,timescale];m$mm1<-m[,timescale]-1
  m[,timescale]<-m[,identifier]<-NULL
  return(m)
}
m<-data_format(e1,"Study_ID","study_month")
m<-data.table(m)
m[, lag1.compl_bw:=c(0, compl_bw[-.N]), by=id]
m[, lag2.compl_bw:=c(0, lag1.compl_bw[-.N]), by=id]
m[, lag1.compl_bw2:=c(0, compl_bw2[-.N]), by=id]
m[, lag2.compl_bw2:=c(0, lag1.compl_bw2[-.N]), by=id]
m[, lag1.compl_bw3:=c(0, compl_bw3[-.N]), by=id]
m[, lag2.compl_bw3:=c(0, lag1.compl_bw3[-.N]), by=id]
m[, lag1.bleed:=c(0, bleed[-.N]), by=id]
m[, lag1.vaginal:=c(0, vaginal[-.N]), by=id]
m[, lag1.gastro:=c(0, gastro[-.N]), by=id]
m[, lag1.conceived:=c(0, conceived[-.N]), by=id]

# look at the data
head(m,5);tail(m,5)

# randomization indicator
R<-as.matrix(as.numeric(m$r));colnames(R)<-c("R")

#baseline variables
V<-as.matrix(cbind(as.factor(m$site),scale(as.numeric(m$rd)),scale(as.numeric(m$age)),
                   as.factor(m$income),as.factor(m$white),as.factor(m$education),as.factor(m$marital),as.factor(m$employed)))
colnames(V)<-c("V1","V2","V3","V4","V5","V6","V7","V8")

#exposure
X<-as.matrix(as.numeric(m$compl_bw3>=.77));colnames(X)<-c("X")
Xl<-as.matrix(as.numeric(m$lag1.compl_bw3>=.77));colnames(Xl)<-c("Xl")

#time-varying confounder 1
Z<-as.matrix(as.numeric(m$conceived>0));colnames(Z)<-c("Z")
Zl<-as.matrix(as.numeric(m$lag1.conceived>0));colnames(Zl)<-c("Zl")

#time-varying confounder 2
B<-as.matrix(as.numeric(m$bleed>0|m$vaginal>0));colnames(B)<-c("B")
Bl<-as.matrix(as.numeric(m$lag1.bleed>0|m$lag1.vaginal>0));colnames(Bl)<-c("Bl")

#time-varying confounder 3
N<-as.matrix(as.numeric(m$gastro>0));colnames(N)<-c("N")
Nl<-as.matrix(as.numeric(m$lag1.gastro>0));colnames(Nl)<-c("Nl")

#censoring and outcomes
m$y<-as.numeric(m$y>0)*m$status;table(m$y)
C<-as.matrix(as.numeric(m$y==4));colnames(C)<-c("C")
S<-as.matrix(as.numeric(m$y==1));colnames(S)<-c("S")
D<-as.matrix(as.numeric(m$y==3));colnames(D)<-c("D")
Y<-as.matrix(as.numeric(m$y==2));colnames(Y)<-c("Y")
j<-as.matrix(as.numeric(ifelse(m$mm>15,15,m$mm)));colnames(j)<-c("j")

#id variable
id<-as.matrix(m$id);colnames(id)<-c("id")

## FIT MODELS FOR LATER PREDICTION
cores<-detectCores()-3

## general function to fit superlearner
# sequential
SL.model <- function(dependent,independent,id,library,cores) {
  output <- SuperLearner(Y=dependent, X=independent,id=id,
                             SL.library=SL.library,
                             method="method.NNLS",family=binomial())
  return(output)
}

## exposure
a<-data.frame(V,R,Xl,Bl,Nl,j)
fitX<-SL.model(dependent=X,independent=a,
               id=id,library=SL.library,
               cores=cores)
fitX

## time-varying confoudner 1
a<-data.frame(V,R,X,Xl,Bl,Nl,j)
fitB<-SL.model(dependent=B,independent=a,
               id=id,library=SL.library,
               cores=cores)
fitB

#time-varying confounder 2
a<-data.frame(V,R,X,Xl,B,Bl,Nl,j)
fitN<-SL.model(dependent=N,independent=a,
               id=id,library=SL.library,
               cores=cores)
fitN

#time-varying confounder 3
a<-data.frame(V,R,X,Xl,B,Bl,N,Nl,Zl,j)[Zl==0,]
a$Zl<-NULL
fitZ<-SL.model(dependent=Z[Zl==0],independent=a,
               id=id[Zl==0],library=SL.library,
               cores=cores)
fitZ

#censoring
a<-data.frame(V,R,X,Xl,B,Bl,N,Nl,Z,j)
fitC<-SL.model(dependent=C,independent=a,
               id=id,library=SL.library,
               cores=cores)
fitC

#outcome 1
a<-data.frame(V,R,X,Xl,B,Bl,N,Nl,Z,j)[Z==0,]
a$Z<-NULL
fitS<-SL.model(dependent=S[Z==0],independent=a,
               id=id[Z==0],library=SL.library,
               cores=cores)
fitS

#outcome 2
a<-data.frame(V,R,X,Xl,B,Bl,N,Nl,Z,j)[Z==1,]
head(a)
a$Z<-NULL
fitD<-SL.model(dependent=D[Z==1],independent=a,
               id=id[Z==1],library=SL.library,
               cores=cores)
fitD

#outcome 3
a<-data.frame(V,R,X,Xl,B,Bl,N,Nl,Z,j)[Z==1,]
a$Z<-NULL
fitY<-SL.model(dependent=Y[Z==1],independent=a,
               id=id[Z==1],library=SL.library,
               cores=cores)
fitY

## PREDICT FOLLOW-UP BASED ON G FORMULA
# sample baseline variables and time-varying variables at first time-point
VV<-data.frame(id,V,R,X,B,N,Z,C,S,D,Y);head(VV,5)
V_first<-VV[!duplicated(VV$id),];head(V_first,5)
index<-strata(V_first,stratanames=c("R"), size=c(montecarlo/2,montecarlo/2), method="srswr");tail(index)
Vpp<-getdata(V_first,index)
drop<-c("ID_unit","Prob","Stratum")
Vpp<-Vpp[,!names(Vpp) %in% drop]
Vpp$id<-1:montecarlo
row.names(Vpp)<-NULL
# initiate dataset
gdat<-NULL
# function for predicting from empirical models
#############
# tested to see if faster, but it is not:
# Unit: milliseconds
# expr      min       lq     mean   median       uq      max neval cld
# ppFunc(fitX, cbind(Vpp[ii, ], Xl = 0, Nl = 0, Bl = 0, j = 2)) 16.32182 20.95740 23.64964 21.72819 22.74474 199.2918  1000   a
# pFunc(fitX, cbind(Vpp[ii, ], Xl = 0, Nl = 0, Bl = 0, j = 2)) 14.75243 21.23699 24.15593 22.04710 23.28962 204.9459  1000   a
#############
# ppFunc <- function (object, newdata){
#   k <- length(object$libraryNames)
#   whichLibrary <- which(object$coef > 0)
#   predY <- matrix(0, nrow = nrow(newdata), ncol = k)
#   for (mm in whichLibrary) {
#     family <- object$family
#     predY[, mm] <- do.call("predict", list(object = object$fitLibrary[[mm]],
#                                           newdata = data.frame(newdata), family = family))
#   }
#   getPred <- object$method$computePred(predY = predY, coef = object$coef,
#                                        control = object$control)
#   out <- list(pred = getPred)
#   return(out)
# }

# prediction function, which takes data and SL object to output predicted value
pFunc<-function(mod,ndat){as.numeric(predict(mod,newdata=ndat,onlySL=T)$pred>runif(1))}
# parametric g formula function
pgf<-function(ii,mc_data,randomization=NULL,exposure=NULL,censoring=NULL){
  # select first observation from baseline data
  Vp<-mc_data[ii,]
  # assign randomization indicator
  Rp<-Vp$R
  # length of follow-up
  length<-15
  # initiate variables to be simulated
  Xp<-Bp<-Np<-Zp<-Cp<-Sp<-Dp<-Yp<-mm<-numeric()
  Xp[1]<-Vp$X;Bp[1]<-Vp$B;Np[1]<-Vp$N;Zp[1]<-Vp$Z;
  Cp[1]<-Vp$C;Sp[1]<-Vp$S;Dp[1]<-Vp$D;Yp[1]<-Vp$Y;
  mm[1]<-j<-1;Vp<-Vp[,names(Vp) %in% c("id","V1","V2","V3","V4","V5","V6","V7","V8")]
  id<-Vp$id;Vp$id<-NULL
  # first time-point is sampled empirical data. for second time-point onwards, simulate
  for(j in 2:length){
    # if no terminal events, then simulate
    if(Cp[j-1]==0&Sp[j-1]==0&Dp[j-1]==0&Yp[j-1]==0) {
      # exposure
      dXp<-data.frame(Vp,R=Rp,Xl=Xp[j-1],Bl=Bp[j-1],Nl=Np[j-1],j)
      Xp[j]<-pFunc(fitX,dXp)
      # time-varying confounder 1
      dBp<-data.frame(Vp,R=Rp,X=Xp[j],Xl=Xp[j-1],Bl=Bp[j-1],Nl=Np[j-1],j)
      Bp[j]<-pFunc(fitB,dBp)
      # time-varying confounder 2
      dNp<-data.frame(Vp,R=Rp,X=Xp[j],Xl=Xp[j-1],B=Bp[j],Bl=Xp[j-1],Nl=Xp[j-1],j)
      Np[j]<-pFunc(fitN,dNp)
      # time-varying confounder 3
      dZp<-data.frame(Vp,R=Rp,X=Xp[j],Xl=Xp[j-1],B=Bp[j],Bl=Bp[j-1],N=Np[j],Nl=Np[j-1],j)
      Zp[j]<-ifelse(Zp[j-1]==0,
                    pFunc(fitZ,dZp),
                    1)
      # outcome 1
      dCp<-data.frame(Vp,R=Rp,X=Xp[j],Xl=Xp[j-1],B=Bp[j],Bl=Bp[j-1],N=Np[j],Nl=Xp[j-1],Z=Zp[j],j)
      Cp[j]<-pFunc(fitC,dCp)
      # outcome 2
      dSp<-data.frame(Vp,R=Rp,X=Xp[j],Xl=Xp[j-1],B=Bp[j],Bl=Xp[j-1],N=Np[j],Nl=Xp[j-1],j)
      if(j<=9){
        Sp[j]<-ifelse(Zp[j]==0&Cp[j]==0&j>4,
                      pFunc(fitS,dSp),
                      0)
      } else{
        Sp[j]<-ifelse(Zp[j]==0&Cp[j]==0,1,0)
      }
      # outcome 3
      dDp<-data.frame(Vp,R=Rp,X=Xp[j],Xl=Xp[j-1],B=Bp[j],Bl=Bp[j-1],N=Np[j],Nl=Np[j-1],j)
      Dp[j]<-ifelse(Cp[j]==0&Sp[j]==0&Zp[j]==1&j>1&j<11,
                    pFunc(fitD,dDp),
                    0)
      # outcome 4
      dYp<-data.frame(Vp,R=Rp,X=Xp[j],Xl=Xp[j-1],B=Bp[j],Bl=Bp[j-1],N=Np[j],Nl=Np[j-1],j)
      Yp[j]<-ifelse(Cp[j]==0&Dp[j]==0&Sp[j]==0&Zp[j]==1&j>7,
                    pFunc(fitY,dYp),
                    0)
    } else {
      break
    }
    # time-point
    mm[j]<-j
  }
  gdat<-as.data.frame(cbind(id,mm,Vp,Rp,Bp,Xp,Zp,Sp,Cp,Dp,Yp))
}

cl <- makeCluster(cores)
registerDoParallel(cl)
clusterCall(cl, function(x) .libPaths(x), .libPaths())
time<-NULL
strt<-Sys.time()
r0<-foreach(i=1:montecarlo,
            .options.RNG=123,
            .packages=c("SuperLearner","glmnet"),
            .export = c("gdat"),
            .errorhandling="pass") %dorng% {
              simp <- pgf(ii=i,mc_data=Vpp)
              simp
            }
stopCluster(cl)
time<-Sys.time()-strt
time

tail(r0)
g_form<-do.call(rbind.data.frame,r0);row.names(r0)<-NULL
head(g_form,10);nrow(g_form)

g0<-ddply(g_form, .(id), function(x) x[nrow(x),]);head(g0)
g0$y<-ifelse(g0$Cp==1,3,ifelse(g0$Sp==1,0,ifelse(g0$Dp==1,2,1)))
table(g0$y)
## STOP
## ## be sure to save both g0 and g_form to file
write.table(g_form, "NatCours1_09Mar17.txt", sep="\t")
write.table(g0, "NatCours2_09Mar17.txt", sep="\t")
