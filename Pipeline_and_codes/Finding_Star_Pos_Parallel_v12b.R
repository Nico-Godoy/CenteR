
memory.limit(16384)

cl <- makeCluster(1)
registerDoParallel(cl)
stopCluster(cl)
cl <- makeCluster(Cluster_number)
registerDoParallel(cl)

work.area<-work.area2

write.table(date(),paste0(work.area,"Date_program_part3.txt"),col.names=FALSE,row.names=FALSE,quote=FALSE,append=FALSE)

###################################################################

A<-read.table(paste0(work.area,"Center_Scs_new2.txt"),header=TRUE)

Anan<-which(is.na(A$Y)==TRUE | is.na(A$X)==TRUE | is.na(A$Yc)==TRUE | is.na(A$Xc)==TRUE)
if(length(Anan)!=0){
  A<-A[-Anan,]
}

B<-read.table(paste0(work.area,"Total_Final_Info_Red.txt"),header=TRUE)
DIMZ<-median(B$DimZ)
# Estimado la distancia entre AGPM y Circle
S<-read.table(paste0(work.area,"Match_Sky_center_new2.txt"),header=TRUE)

Snan<-which(is.na(S$Y)==TRUE | is.na(S$X)==TRUE | is.na(S$Yc)==TRUE | is.na(S$Xc)==TRUE)

if(length(Snan)!=0){
  S<-S[-Snan,]
}


# X: AGPM Xc:Circle
Xt<- S$Xc-S$X
DXt<-sqrt( (S$Err_Xc)**2 + (S$Err_X)**2)
Yt<- S$Yc-S$Y
DYt<-sqrt( (S$Err_Yc)**2 + (S$Err_Y)**2)

# Estadisticamente

XT<-median(Xt) # -10.45105
DXT<-0.67499*mad(Xt)  #   0.02003
YT<-median(Yt) # +11.71069
DYT<-0.67499*mad(Yt)  #   0.01092

write.table(matrix(c(XT,DXT,YT,DYT),nr=1),paste0(work.area,"Transf.txt"),col.names=FALSE,row.names=FALSE)

## The object to test:

SupM<-matrix(nr=0,nc=dim(A)[2]+dim(B)[2]+19+2)

Bname<-as.character(unique(B$NewName))

write.table(paste0("_",0,"_",length(Bname)),"Where.txt",col.names=FALSE,row.names=FALSE,append=FALSE)

VALS<-foreach(IU=1:length(Bname), .combine=rbind,.packages=c('foreach','mvtnorm','FITSio','fields','astro')) %dopar%  {
#VALS<-foreach(IU=129:130, .combine=rbind,.packages=c('foreach','mvtnorm','FITSio','fields','astro')) %dopar% {

nam<-which(B$NewName==Bname[IU])
nam2<-which(A$Name==paste0("Scs_",as.character(B$FILE[nam[1]])))

B0<-B[nam,]
if(length(nam2)!=0){

A0<-A[nam2,]

if(A0$Xc>2024){
  if(nam2>2 & nam2<(dim(A)[1]-1)){
    A0$Xc<-mean(c(A$Xc[(nam2-2):(nam2-1)],A$Xc[(nam2+1):(nam2+2)]),na.rm=TRUE)
  }else if(nam2<2){
    A0$Xc<-mean(A$Xc[(nam2+(3:4))],na.rm=TRUE)
  }else{
    A0$Xc<-mean(A$Xc[(nam2-(3:4))],na.rm=TRUE)
  }
}

if(A0$Yc>2024){
  if(nam2>2 & nam2<(dim(A)[1]-1)){
    A0$Yc<-mean(c(A$Yc[(nam2-2):(nam2-1)],A$Yc[(nam2+1):(nam2+2)]),na.rm=TRUE)
  }else if(nam2<2){
    A0$Yc<-mean(A$Yc[(nam2+(3:4))],na.rm=TRUE)
  }else{
    A0$Yc<-mean(A$Yc[(nam2-(3:4))],na.rm=TRUE)
  }
}


if(is.na(A0$Xc)==FALSE){

if(length(B0$DimZ)>50 ){
DIMZ<-B0$DimZ[1]
B0<-B0[1:DIMZ,]


################################################################


Fun<-function(p,X,Y,Z,Xc,Yc,TF){
 if(p[1]>min(X) & p[1]<max(X) & p[2]>min(Y) & p[2]<max(Y) & sqrt(p[3]**2 + p[4]**2)>p[6] &
    p[3]<100*length(X)/2 &  p[4]<100*length(Y)/2 & p[7]>0 & p[8]>=0 &  abs(p[5])<=360 &
     p[6]<=2.15*100 & p[6]>=0.55*100 & p[8]/p[7]<0.999){# &
## @ 06-08-2019 changed the limit for p[6] from 0.5 to 1.5
   xy<-expand.grid(X,Y)
   MeanStar<-c(p[1],p[2])
   SigStar0<-matrix(c( (p[3]**2)/1e4 ,0,0,(p[4]**2)/1e4),nc=2)
##
   theta<-(p[5]*pi/180)
   RotM<- matrix(c(cos(theta),sin(theta),-sin(theta),cos(theta)),nrow=2)
   SigStar<-RotM %*% SigStar0 %*% t(RotM)
##
   MeanCor<-c(Xc[1],Yc[1])
   SigCor<-matrix(c(p[6]**2,0,0,p[6]**2),nc=2)/1e4
   VStar<-matrix(dmvnorm(xy,MeanStar,SigStar),nc=length(X),nr=length(Y))
   VCor<-matrix(dmvnorm(xy,MeanCor,SigCor),nc=length(X),nr=length(Y))
   Z<-Z/abs(sum(Z))
   V<-p[7]*VStar - p[8]*VCor + median(Z[TF==TRUE])
   S<-sum(( (Z-V)**2),na.rm=T)*100
  return(S)
 }else{
  return(1e298)
 }
}

fun<-function(p,X,Y,Z,Xc,Yc,TF){
   xy<-expand.grid(X,Y)
   MeanStar<-c(p[1],p[2])
   SigStar0<-matrix(c( (p[3]**2)/1e4 ,0,0,(p[4]**2)/1e4),nc=2)
##
   theta<-p[5]*pi/180
   RotM<- matrix(c(cos(theta),sin(theta),-sin(theta),cos(theta)),nrow=2)
   SigStar<-RotM %*% SigStar0 %*% t(RotM)
##
   MeanCor<-c(Xc[1],Yc[1])
   SigCor<-matrix(c(p[6]**2,0,0,p[6]**2),nc=2)/1e4
   VStar<-matrix(dmvnorm(xy,MeanStar,SigStar),nc=length(X),nr=length(Y))
   VCor<-matrix(dmvnorm(xy,MeanCor,SigCor),nc=length(X),nr=length(Y))
   Z<-Z/abs(sum(Z))
   V<-p[7]*VStar - p[8]*VCor + median(Z[TF==TRUE])
   return(V)
}

Xc<-A0$Xc - XT
Yc<-A0$Yc - YT
DXc<-sqrt( DXT**2 + A0$Err_Xc**2  )
DYc<-sqrt( DYT**2 + A0$Err_Yc**2  )
R<-11
Name0<-strsplit(paste0("/Scs_",as.character(B$FILE[nam[1]])),".fits")[[1]]
Img<-readFITS(paste0(work.area,"Data/",Name0,"_DF.fits"))$imDat

Z0<-readFITS(paste0(work.area,"Data/",as.character(A0$Name)))$imDat[floor(Xc-R):floor(Xc+R),floor(Yc-R):floor(Yc+R)]
Z<-Img[floor(Xc-R):floor(Xc+R),floor(Yc-R):floor(Yc+R),]

X0<-Xc-R
Y0<-Yc-R

X<-floor(Xc-R):floor(Xc+R)
Y<-floor(Yc-R):floor(Yc+R)

if(min(X)<0){ write.table(IU,paste0(work.area,"/meme.txt")) }


TF<-matrix(TRUE,nr=length(X),nc=length(Y))
xyF<-expand.grid(1:length(X),1:length(Y))
tf<-which(sqrt( (xyF[,1]-(Xc-X[1]+1))**2 + (xyF[,2]-(Yc-Y[1]+1))**2)<=5  )
TF[tf]<-FALSE

####

 ANG0<-seq(0,360,len=10)
 posX<-seq(Xc-1,Xc+1,len=10)
 posY<-seq(Yc-1,Yc+1,len=10)
 SgXY<-seq(2.8,4.8,len=10)
 Sg0<-seq(0.61,2.0,len=10)
 Int0<-seq(0.3e1,5e1,len=10)
 Int0c<-Int0*0.35

 XYTpos<-as.matrix(expand.grid(posX,posY))
 Nlen<-15
 aq12<- XYTpos[sample( 1:length(XYTpos[,1]),Nlen,rep=TRUE),]
 aq34<- matrix(sample(SgXY,Nlen*2,rep=TRUE),nc=2)
 aq5<- matrix(sample(ANG0,Nlen,rep=TRUE),nc=1)
 aq6<- matrix(sample(Sg0,Nlen,rep=TRUE) ,nc=1)
 aq7<- matrix(sample(Int0,Nlen,rep=TRUE),nc=1)
 aq8<- aq7*0.4

 XYTpos<-cbind(aq12,aq34,aq5,aq6,aq7,aq8)
 
 XYTpos<- matrix(as.vector(XYTpos),nc=8)

####

P0<-matrix(nr=0,nc=8)
P0<-rbind(P0,c( Xc,Yc, 3.1,3.1 , 0, 1.8, 1e1,0.8e1))
P0<-rbind(P0,c( Xc,Yc, 3.1,3.1 , 0, 1.8, 1e1,0.8e1))

P0<-rbind(P0,c( Xc+.5,Yc+.5, 3.1,5.1 , 30, 0.6, 1e1,0.3e1))
#P0<-rbind(P0,c( Xc+.5,Yc-.5, 4.1 , -5e-4, 1.5, 1e1,0.3e1))
P0<-rbind(P0,c( Xc-.5,Yc+-.5, 4.2,5.1 , 330,0.6, 1e1,0.3e1))
#P0<-rbind(P0,c( Xc-.5,Yc-.5, 4.3 , -5e-4, 1.5, 1e1,0.3e1))

P0<-rbind(P0,c( Xc+1,Yc-1, 3.1,5.1 , 60, 1.9, 1e1,0.3e1))
#P0<-rbind(P0,c( Xc+1,Yc+1.1, 3.1 , 5e-3, 1.9, 1e1,0.3e1))
P0<-rbind(P0,c( Xc+1,Yc+1, 2.2,3.2 , 300, 1.9, 1e1,0.3e1))
#P0<-rbind(P0,c( Xc+1,Yc-1, 3.3 , 5e-3, 1.9, 1e1,0.3e1))

P0<-rbind(P0,c( Xc+3.1,Yc, 4.4 ,4.4 , 90, 2, 2e1,0.3e1))
#P0<-rbind(P0,c( Xc-3.2,Yc, 4.6 , 0.001, 1, 2e1,0.3e1))
P0<-rbind(P0,c( Xc,Yc+3.1, 3.5,3.5  , 120, 2, 2e1,0.3e1))
#P0<-rbind(P0,c( Xc,Yc-3.2, 4.5 , 0.001, 1, 2e1,0.3e1))

P0<-rbind(P0,c( Xc+3.1,Yc+3.2, 3.7,3.7 , 240, 0.9, 0.7e1,0.3e1))
#P0<-rbind(P0,c( Xc+3.2,Yc+3.1, 4 , 0.001, 1, 0.7e1,0.3e1))
P0<-rbind(P0,c( Xc+3.2,Yc+3.1,4.7 ,4.7, 180, 0.9, 0.7e1,0.3e1))
#P0<-rbind(P0,c( Xc-3.1,Yc-3.2, 5.5 , 0.001, 1, 0.7e1,0.3e1))

P0<-rbind(P0,c( Xc,Yc-6.1 ,4.7 ,4.7 , 45, 1.2, 1e1,0.6e1))
#P0<-rbind(P0,c( Xc,Yc+6.1, 5 , 0.00001, 1.8, 1e1,0.6e1))
#P0<-rbind(P0,c( Xc-6.2,Yc, 6.1 , 0.00001, 1.8, 1e1,0.6e1))
P0<-rbind(P0,c( Xc+6.2,Yc, 4.3 , 4.3,315, 1.2, 1e1,0.6e1))

P0<-rbind(P0,c( Xc-4,Yc-4.1, 4.1, 4.1, 25, 1.7, 1e1,0.6e1))
P0<-rbind(P0,c( Xc-4,Yc-4.1, 4.1,4.1 , 15, 1.7, 1e1,0.6e1))
#P0<-rbind(P0,c( Xc-6,Yc+6.1, 6  , -0.1e-3, 1.7, 1e1,0.6e1))
#P0<-rbind(P0,c( Xc+6.2,Yc-6, 6.2 , -0.1e-3, 1.7, 1e1,0.6e1))
P0<-rbind(P0,c( Xc+.2,Yc-.7, 3.8 ,3.8, 115, 0.7, 1e1,0.6e1))
P0<-rbind(P0,c( Xc-.7,Yc+.2, 3.8 ,3.8, 105, 0.7, 1e1,0.6e1))

#P0<-P0[c(1,2,3,4,5,10,11,12,13,18,19,20,21,6,7,8,9,14,15,16,17,22,23,24,25),]
#P0[,c(3,4,5,6)]<-P0[,c(3,4,5,6)]*100

###
P0<-XYTpos
###

P0[,c(3,4,5,6)]<-P0[,c(3,4,5,6)]*100
P0[,5]<-P0[,5]/100

Tol<-seq(1e-6,0.1e-3,len=dim(P0)[1])
STM<-seq(10000,300000,len=dim(P0)[1])
TolG<-seq(1e-7,1e-3,len=dim(P0)[1])

Nind<-8

SMI<-matrix(nr=DIMZ,nc=25)
for(i in 1:DIMZ){
  N<-i
  SS<-1
  XM<-which(Z[,,N]==max(Z[,,N]),arr.ind=TRUE)
  XM<-c(X[XM[1]],Y[XM[2]])
  P1<-P0[SS,]
  P1[1:2]<-XM[1:2]
  P1[3:4]<-P1[3:4]/2

  MIN<-0

#  while(SS<(dim(P0)[1]+1) ){
  for(SS in 1:(dim(P0)[1]+1)){
   if(SS<=dim(P0)[1]){
    if(class(try(nlm(Fun,p=P0[SS,],X=X,Y=Y,Z=Z[,,N],Xc=c(Xc,DXc),Yc=c(Yc,DYc),TF=TF,hessian=TRUE,iterlim = 8000,
         steptol = Tol[SS],gradtol = TolG[SS],stepmax=STM[SS])))!='try-error'){
      res<-nlm(Fun,p=P0[SS,],X=X,Y=Y,Z=Z[,,N],Xc=c(Xc,DXc),Yc=c(Yc,DYc),TF=TF,hessian=TRUE,iterlim = 8000,
          steptol = Tol[SS],gradtol = TolG[SS],stepmax=STM[SS])
      #SS<-20000000
    }else{
      res<-list(est=c(9999+(1:Nind)*0),heSSiian=matrix(0,nr=Nind,nc=Nind),minimum=Inf)
      cat("Not found a non-singularity solution ...\n")
    }
      #SS<-SS+1}
   }else{
    if(class(try(nlm(Fun,p=res$est,X=X,Y=Y,Z=Z[,,N],Xc=c(Xc,DXc),Yc=c(Yc,DYc),TF=TF,hessian=TRUE,iterlim = 8000,
         steptol = 1e-5,gradtol = 1e-5,stepmax=1000)))!='try-error'){
      res<-nlm(Fun,p=res$est,X=X,Y=Y,Z=Z[,,N],Xc=c(Xc,DXc),Yc=c(Yc,DYc),TF=TF,hessian=TRUE,iterlim = 8000,
         steptol = 1e-5,gradtol = 1e-5,stepmax=1000)
      #SS<-20000000
    }else if(class(try(nlm(Fun,p=res$est,X=X,Y=Y,Z=Z[,,N],Xc=c(Xc,DXc),Yc=c(Yc,DYc),TF=TF,hessian=TRUE,iterlim = 8000,
         steptol = 1e-3,gradtol = 1e-3,stepmax=1000)))!='try-error' & i!=1){
      res<-nlm(Fun,p=res$est,X=X,Y=Y,Z=Z[,,N],Xc=c(Xc,DXc),Yc=c(Yc,DYc),TF=TF,hessian=TRUE,iterlim = 8000,
         steptol = 1e-3,gradtol = 1e-3,stepmax=1000)
      #SS<-20000000
    }else{
      res<-list(est=c(9999+(1:Nind)*0),heSSiian=matrix(0,nr=Nind,nc=Nind),minimum=Inf)
      cat("Not found a non-singularity solution ...\n")
      #SS<-20000000
    }
   }
   MIN[SS]<-res$min
  }

  minV<-which(MIN==min(MIN))

  SSi<-minV[1]
  
   if(SSi<=dim(P0)[1]){
    if(class(try(nlm(Fun,p=P0[SSi,],X=X,Y=Y,Z=Z[,,N],Xc=c(Xc,DXc),Yc=c(Yc,DYc),TF=TF,hessian=TRUE,iterlim = 8000,
         steptol = Tol[SSi],gradtol = TolG[SSi],stepmax=STM[SSi])))!='try-error'){
      res<-nlm(Fun,p=P0[SSi,],X=X,Y=Y,Z=Z[,,N],Xc=c(Xc,DXc),Yc=c(Yc,DYc),TF=TF,hessian=TRUE,iterlim = 8000,
          steptol = Tol[SSi],gradtol = TolG[SSi],stepmax=STM[SSi])
      SSi<-20000000
    }else{
      SSi<-SSi+1} 
   }else{
    if(class(try(nlm(Fun,p=res$est,X=X,Y=Y,Z=Z[,,N],Xc=c(Xc,DXc),Yc=c(Yc,DYc),TF=TF,hessian=TRUE,iterlim = 8000,
         steptol = 1e-5,gradtol = 1e-5,stepmax=1000)))!='try-error'){
      res<-nlm(Fun,p=res$est,X=X,Y=Y,Z=Z[,,N],Xc=c(Xc,DXc),Yc=c(Yc,DYc),TF=TF,hessian=TRUE,iterlim = 8000,
         steptol = 1e-5,gradtol = 1e-5,stepmax=1000)
      SSi<-20000000
    }else if(class(try(nlm(Fun,p=res$est,X=X,Y=Y,Z=Z[,,N],Xc=c(Xc,DXc),Yc=c(Yc,DYc),TF=TF,hessian=TRUE,iterlim = 8000,
         steptol = 1e-3,gradtol = 1e-3,stepmax=1000)))!='try-error' & i!=1){
      res<-nlm(Fun,p=res$est,X=X,Y=Y,Z=Z[,,N],Xc=c(Xc,DXc),Yc=c(Yc,DYc),TF=TF,hessian=TRUE,iterlim = 8000,
         steptol = 1e-3,gradtol = 1e-3,stepmax=1000)
      SSi<-20000000
    }else{
      res<-list(est=c(9999+(1:Nind)*0),heSSiian=matrix(0,nr=Nind,nc=Nind),minimum=Inf)
      cat("Not found a non-singularity solution ...\n")
      SSi<-20000000
    }
   }


  if(class(try(solve(res$hessian)))=='matrix'){
     Err<-sqrt(abs(diag(Nind*(res$minimum)/(length(X)*length(Y) - Nind)*solve(res$hessian))))*1.96
  }else if(class(try(solve(res$hessian,tol=1e-8)))=='matrix'){
     Err<-sqrt(abs(diag(Nind*(res$minimum)/(length(X)*length(Y) - Nind)*solve(res$hessian,tol=1e-8))))*1.96
  }else if(class(try(solve(res$hessian,tol=1e-4)))=='matrix'){
     Err<-sqrt(abs(diag(Nind*(res$minimum)/(length(X)*length(Y) - Nind)*solve(res$hessian,tol=1e-4))))*1.96
  }else if(class(try(solve(res$hessian,tol=1e-2)))=='matrix'){
     Err<-sqrt(abs(diag(Nind*(res$minimum)/(length(X)*length(Y) - Nind)*solve(res$hessian,tol=1e-1))))*1.96
  }else{
     Err<-0*(1:Nind)
  }

  Zkk<-Z[,,N]
  bkgK<-median(Zkk[TF==TRUE])
SMI[i,]<-c(i,res$est[1],Err[1],res$est[2],Err[2],res$est[3],Err[3],res$est[4],Err[4],res$est[5],Err[5],
           res$est[6],Err[6],bkgK,0,res$est[7]*abs(sum(Z[,,N])),Err[7]*abs(sum(Z[,,N])),res$est[8]*abs(sum(Z[,,N])),Err[8]*abs(sum(Z[,,N])),
            Xc[1],DXc[1], Yc[1],DYc[1],res$minimum,sd(Zkk[TF==TRUE]))

cat(paste0("--- ",i,"/ ",dim(SMI)[1],"  ---",IU,"/",length(Bname),"\n"))

}

Dist<-sqrt( (Xc-SMI[,2])**2 + (Yc-SMI[,4])**2   )
#E_Dist<- (1/Dist)*sqrt( (DXc**2 + SMI[,3]**2)*(Xc-SMI[,2])**2 + (DYc**2 + SMI[,5]**2)*(Yc-SMI[,4])**2   ) 

ASg<-sqrt( (SMI[,6]**2 + SMI[,8]**2)/2  )
#E_ASg<- (1/(2*ASg))*sqrt(  (SMI[,6]*SMI[,7])**2 +  (SMI[,8]*SMI[,9])**2)

R<-sqrt( SMI[,6] + SMI[,10]*sqrt(SMI[,6]*SMI[,8]))/sqrt( SMI[,6] - SMI[,10]*sqrt(SMI[,6]*SMI[,8]) )

cat("--is ok -TP...\n")
#SPF<-cbind(B0,SMI,matrix(c(Dist,E_Dist,ASg,E_ASg,R),nc=5),t(matrix(A0,nc=dim(SMI)[1],nr=dim(A0)[2]))   )
SPF<-cbind(B0,SMI,matrix(c(Dist,ASg,R),nc=3),t(matrix(A0,nc=dim(SMI)[1],nr=dim(A0)[2]))   )

cat("--is ok...\n")
write.table(paste0("_",IU,"_",length(Bname)),"Where.txt",col.names=FALSE,row.names=FALSE,append=TRUE)
 return(SPF)
}else{ # if of dim-Z
   kk0<-cbind(B0[1,],matrix(99999,nr=1,nc=25),matrix(c(99999,99999,99999),nc=3),t(matrix(99999,nc=1,nr=dim(A)[2])))
   SPF<-matrix(nr=0,nc=length(kk0[1,]))
  for(ITU in 1:10){
   SPF<-rbind(SPF,kk0,kk0,kk0,kk0,kk0,kk0,kk0,kk0,kk0,kk0)
  }

 write.table(paste0("_",IU,"_",length(Bname)),"Where.txt",col.names=FALSE,row.names=FALSE,append=TRUE)
 return(SPF)
}
 #   p1<-c( Xc+2,Yc+0.5 , 7 , 7 ,0,1.5,1.5,1e1,0.8e1)
}else{  # if of A0$Xc
   kk0<-cbind(B0[1,],matrix(99999,nr=1,nc=25),matrix(c(99999,99999,99999),nc=3),t(matrix(99999,nc=1,nr=dim(A)[2])))
   SPF0<-matrix(nr=0,nc=length(kk0[1,]))
  for(ITU2 in 1:10){
   SPF0<-rbind(SPF0,kk0,kk0,kk0,kk0,kk0,kk0,kk0,kk0,kk0,kk0)
  }
   return(SPF0)
}

}else{  # if of length(Nam2)
   kk0<-cbind(B0[1,],matrix(99999,nr=1,nc=25),matrix(c(99999,99999,99999),nc=3),t(matrix(99999,nc=1,nr=dim(A)[2])))
   SPF0<-matrix(nr=0,nc=length(kk0[1,]))
  for(ITU2 in 1:10){
   SPF0<-rbind(SPF0,kk0,kk0,kk0,kk0,kk0,kk0,kk0,kk0,kk0,kk0)
  }
   return(SPF0)
}
} # foreach of IU

cat("All is ok at this point\n")
headers<-c(names(B), "index","Sx","E_Sx","Sy","E_Sy","SDx","E_SDx","SDy","E_SDy","Ang","E_Ang","CoS",
        "E_CoS","Bkg","E_Bkg","Nor1","E_Nor1","Nor2","E_Nor2","CorX","E_CorX","CorY","E_CorY","RES","BkgSD","DisCS","EffSD","Ratio",names(A) )
colnames(VALS)<-headers

VALS2<-as.matrix(VALS)

write.table(VALS2,paste0(work.area,"/Final_table_complete_information_v7.txt"),col.names=TRUE,row.names=FALSE,quote=FALSE)
stopCluster(cl) # ***
#write.table(VALS,paste0(work.area,"/Final_table_complete_information.txt"),col.names=TRUE,row.names=FALSE,quote=FALSE)


write.table(date(),paste0(work.area,"Date_program_part3.txt"),col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE)




# bsub -q server_name -R "rusage[mem=requested_memory]" "Rscript script_name.R"





