#################
source(file="~/Dropbox/Documents/R/ipak.R")
packages <- c("SuperLearner","foreach","doParallel","doRNG","parallel","beepr","data.table","plyr","cmprsk","VIM","mice",
              "lattice","xtable","VIM","sampling","betareg","boot","arm","ipred")
ipak(packages)
e1<-read.table("~/Dropbox/Documents/Research/Papers/EAGeR_GFormula/imputed_eager_monthly_final_03Jan17.txt", sep="\t",header=T)
source(file="~/Dropbox/Documents/Research/Papers/EAGeR_GFormula/EAGeRWrapper.R")
SL.library<-c("SL.mean",#"SL.svm",
              "SL.glm",#"SL.nnet",
              "SL.gam.1","SL.gam.2","SL.gam.3","SL.gam.4",
              #"SL.gbm","SL.gbm.1","SL.gbm.2","SL.gbm.3","SL.gbm.4",
              "SL.glmnet","SL.glmnet.0.25","SL.glmnet.0.5","SL.glmnet.0.75")
SL.library<-c("SL.mean","SL.bayesglm","SL.ipredbagg")
#################

## ARRANGE DATASET & CREATE VARIABLES
m<-e1 
mc_samp<-1000;nrow(m)
m$id<-m$Study_ID;m$Study_ID<-NULL
m$mm<-m$study_month;m$mm1<-m$study_month-1
frac<-sum(as.numeric(m$mm==1))/sum(as.numeric(e1$study_month==1))
m<-data.table(m);head(m)
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

#randomization indicator
R<-as.matrix(as.numeric(m$r));colnames(R)<-c("R")
#baseline variables
V<-as.matrix(cbind(as.factor(m$site),scale(as.numeric(m$rd)),scale(as.numeric(m$age)),
                   as.factor(m$income),as.factor(m$white),as.factor(m$education),as.factor(m$marital),as.factor(m$employed)))
colnames(V)<-c("V1","V2","V3","V4","V5","V6","V7","V8")
#exposure
X<-as.matrix(as.numeric(m$compl_bw3>=1.77));colnames(X)<-c("X")
Xl<-as.matrix(as.numeric(m$lag1.compl_bw3>=1.77));colnames(Xl)<-c("Xl")
#time-varying confounder 1
Z<-as.matrix(as.numeric(m$conceived>1));colnames(Z)<-c("Z")
Zl<-as.matrix(as.numeric(m$lag1.conceived>1));colnames(Zl)<-c("Zl")
#time-varying confounder 2
B<-as.matrix(as.numeric(m$bleed>1|m$vaginal>1));colnames(B)<-c("B")
Bl<-as.matrix(as.numeric(m$lag1.bleed>1|m$lag1.vaginal>1));colnames(Bl)<-c("Bl")
#time-varying confounder 3
N<-as.matrix(as.numeric(m$gastro>1));colnames(N)<-c("N")
Nl<-as.matrix(as.numeric(m$lag1.gastro>1));colnames(Nl)<-c("Nl")
#outcomes
m$y<-as.numeric(m$y>0)*m$status;table(m$y)
C<-as.matrix(as.numeric(m$y==4));colnames(C)<-c("C")
S<-as.matrix(as.numeric(m$y==1));colnames(S)<-c("S")
D<-as.matrix(as.numeric(m$y==3));colnames(D)<-c("D")
Y<-as.matrix(as.numeric(m$y==2));colnames(Y)<-c("Y")
j<-as.matrix(as.numeric(ifelse(m$mm>15,15,m$mm)));colnames(j)<-c("j")
#id variable
id<-as.matrix(m$id);colnames(id)<-c("id")

## FIT MODELS TO DATA FOR LATER PREDICTION 
cores<-detectCores()-1
#exposure
design.fitX<-data.frame(cbind(V,Xl,Bl,Nl,j))
names(design.fitX)<-c("V1","V2","V3","V4","V5","V6","V7","V8","Xl","Bl","Nl","j")
Xm<-as.numeric(X)
idm<-id
cl<-makeCluster(cores,type="FORK")
clusterSetRNGStream(cl, iseed = 123)
fitX <- snowSuperLearner(cluster=cl,Y=Xm, X=design.fitX,id=idm,
                         SL.library=SL.library,
                         method="method.NNLS",family=binomial())
stopCluster(cl)
fitX

## time-varying confoudner 1
design.fitB<-data.frame(cbind(V,X,Xl,Bl,Nl,j))
names(design.fitB)<-c("V1","V2","V3","V4","V5","V6","V7","V8","X","Xl","Bl","Nl","j")
Bm<-as.numeric(B)
idm<-id
cl<-makeCluster(cores,type="FORK")
clusterSetRNGStream(cl, iseed = 123)
fitB <- snowSuperLearner(cluster=cl, Y=Bm, X=design.fitB,id=idm,
                         SL.library=SL.library,
                         method="method.NNLS",family=binomial())
stopCluster(cl)
fitB

#time-varying confounder 2
design.fitN<-data.frame(cbind(V,R,X,Xl,B,Bl,Nl,j))
names(design.fitN)<-c("V1","V2","V3","V4","V5","V6","V7","V8","R","X","Xl","B","Bl","Nl","j")
Nm<-as.numeric(N)
idm<-id
cl<-makeCluster(cores,type="FORK")
clusterSetRNGStream(cl, iseed = 123)
fitN <- snowSuperLearner(cluster=cl, Y=Nm, X=design.fitN,id=idm,
                         SL.library=SL.library,
                         method="method.NNLS",family=binomial())
stopCluster(cl)
fitN

#time-varying confounder 3
design.fitZ<-data.frame(cbind(V,R,X,Xl,B,Bl,N,Nl,Zl,j))[Zl==0,]
design.fitZ$Zl<-NULL
names(design.fitZ)<-c("V1","V2","V3","V4","V5","V6","V7","V8","R","X","Xl","B","Bl","N","Nl","j")
Zm<-as.numeric(Z)[Zl==0]
idm<-id[Zl==0]
cl<-makeCluster(cores,type="FORK")
clusterSetRNGStream(cl, iseed = 123)
fitZ <- snowSuperLearner(cluster=cl,Y=Zm, X=design.fitZ,id=idm,
                         SL.library=SL.library,
                         method="method.NNLS",family=binomial())
stopCluster(cl)
fitZ

#outcome 1
design.fitC<-data.frame(cbind(V,R,X,Xl,B,Bl,N,Nl,Z,j))
names(design.fitC)<-c("V1","V2","V3","V4","V5","V6","V7","V8","R","X","Xl","B","Bl","N","Nl","Z","j")
Cm<-C
idm<-id
cl<-makeCluster(cores,type="FORK")
clusterSetRNGStream(cl, iseed = 123)
fitC <- snowSuperLearner(cluster=cl,Y=Cm, X=design.fitC,id=idm,
                         SL.library=SL.library,
                         method="method.NNLS",family=binomial())
stopCluster(cl)
fitC

#outcome 2
design.fitS<-data.frame(cbind(V,R,X,Xl,B,Bl,N,Nl,Z,j))[Z==0,]
design.fitS$Z<-NULL
names(design.fitS)<-c("V1","V2","V3","V4","V5","V6","V7","V8","R","X","Xl","B","Bl","N","Nl","j")
Sm<-as.numeric(S)[Z==0]
idm<-as.numeric(id)[Z==0]
# check subsetting (TRUE x 2)
(nrow(design.fitS)==length(Sm))
(nrow(design.fitS)==length(idm))
cl<-makeCluster(cores,type="FORK")
clusterSetRNGStream(cl, iseed = 123)
fitS <- snowSuperLearner(cluster=cl,Y=Sm, X=design.fitS,id=idm,
                         SL.library=SL.library,
                         method="method.NNLS",family=binomial())
stopCluster(cl)
fitS

#outcome 3
design.fitD<-data.frame(cbind(V,R,X,Xl,B,Bl,N,Nl,Z,j))
names(design.fitD)<-c("V1","V2","V3","V4","V5","V6","V7","V8","R","X","Xl","B","Bl","N","Nl","Z","j")
design.fitD<-subset(design.fitD,Z==1)
design.fitD$Z<-NULL
Dm<-as.numeric(D)[Z==1]
idm<-as.numeric(id)[Z==1]
cl<-makeCluster(cores,type="FORK")
clusterSetRNGStream(cl, iseed = 123)
fitD <- snowSuperLearner(cluster=cl,Y=Dm, X=design.fitD,id=idm,
                          SL.library=SL.library,
                          method="method.NNLS",family=binomial())
stopCluster(cl)
fitD

#outcome 4
design.fitY<-data.frame(cbind(V,R,X,Xl,B,Bl,N,Nl,Z,j))[Z==1,]
design.fitY$Z<-NULL
names(design.fitY)<-c("V1","V2","V3","V4","V5","V6","V7","V8","R","X","Xl","B","Bl","N","Nl","j")
Ym<-as.numeric(Y)[Z==1]
idm<-as.numeric(id)[Z==1]
cl<-makeCluster(cores,type="FORK")
clusterSetRNGStream(cl, iseed = 123)
fitY <- snowSuperLearner(cluster=cl,Y=Ym, X=design.fitY,id=idm,
                         SL.library=SL.library,
                         method="method.NNLS",family=binomial())
stopCluster(cl)
fitY

## PREDICT FOLLOW-UP BASED ON G FORMULA
# montecarlo sample size
montecarlo<-mc_samp
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
write.table(g_form, "g_form_e0_NatCours.txt", sep="\t")
write.table(g0, "g0_e0_NatCours.txt", sep="\t")
