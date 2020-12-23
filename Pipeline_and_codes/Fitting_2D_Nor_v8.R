
work.area<-work.area2
path.data<-paste0(work.area,"Data/")

setwd(path.data)
names<-system("ls -d *_NACO*.fits | grep -v \'DF\'",intern =TRUE)
setwd(work.area)


Fun<-function(p,X,Y,Z){
 if(p[1]>min(X) & p[1]<max(X) & p[2]>min(Y) & p[2]<max(Y) & p[3]>0 & p[4]>0 & p[3]<30 & p[4]<30){
   xy<-expand.grid(X,Y)
   Mean<-c(p[1],p[2])
   Sig<-matrix(c(p[3],0,0,p[4]),nc=2)
   V<-t(matrix(dmvnorm(xy,Mean,Sig),nc=length(X),nr=length(Y)))
   V<-V/sum(V)
   S<-sum(( (Z-V)**2))*1000
  return(S)
 }else{
  return(1e308)
 }
}


Fun2<-function(p,X,Z){
 #print(p)
  if(p[2]>0 & p[4]>0 & p[3]>0 & p[5]>0){
   y<-p[3]*exp( -((X-p[1])**2)/(2*p[2])  ) + p[6]/2
   y2<-p[5]*exp( -((X-p[1])**2)/(2*p[4])  ) + p[6]/2
   S<-sum((Z-y-y2)**2)
   return(S*1e10)
  }else{
   return(1e270)
  }
}


 CENT2<-function(A){
    if(is.na(dim(A)[3])==FALSE){
     E<-A[,,1]
     for(i in 1:length(A[1,,1])){
      for(j in 1:length(A[,1,1])){
       E[j,i]<-sum(A[j,i,])
     }}
    }else{
     E<-A 
    }
    XYap2<-function(E,x,y,fact){
     V1<-0
     alp<-seq(0,pi,len=length(x)/fact)
     x2<-unique(round((mean(x)-1)*(1 + cos(alp))+1))
     y2<-unique(round((mean(y)-1)*(1 + cos(alp))+1))

     NonX<-seq(4,dim(E)[1],by=8)   # (3 + 1) + 7 ...
     NonY<-seq(4,dim(E)[2],by=8)  

  E[NonX,]<- -99999
  E[,NonY]<- -99999
    V1<-0
    V2<-0
    NoU<-which(is.na(E)==TRUE | E==-99999)
    if(length(NoU)!=0){
       Med<-median(E[-NoU],na.rm=T)
       Mad<-mad(E[-NoU],na.rm=T)
       E[NoU]<- -99999
    }else{
       Med<-median(E,na.rm=T)
       Mad<-mad(E,na.rm=T)
    }
   if(dim(A)[1]<900 | dim(A)[2]<900){
    for(i in 1:max(c(length(x2),length(y2)))){
      if(i<=length(x2)){
        V1[i]<-length(which( E[x2[i],]>=max(c(0,Med-10*Mad)) &  E[x2[i],]<=Med+30*Mad  ))
      }
      if(i<=length(y2)){
        V2[i]<-length(which( E[,y2[i]]>=max(c(0,Med-10*Mad)) &  E[,y2[i]]<=Med+30*Mad  ))
      }
    }
   }else{
    for(i in 1:max(c(length(x2),length(y2)))){
      if(i<=length(x2)){
        V1[i]<-length(which( E[x2[i],]>=max(c(0,Med+3*Mad)) &  E[x2[i],]<=Med+30*Mad  ))
      }
      if(i<=length(y2)){
        V2[i]<-length(which( E[,y2[i]]>=max(c(0,Med+3*Mad)) &  E[,y2[i]]<=Med+30*Mad  ))
      }
    }
   }
#    for(i in 1:max(c(length(x2),length(y2)))){
#      if(i<=length(x2)){
#        V1[i]<-length(which( E[x2[i],]>=max(c(0,Med+3*Mad)) &  E[x2[i],]<=Med+30*Mad  ))
#      }
#      if(i<=length(y2)){
#        V2[i]<-length(which( E[,y2[i]]>=max(c(0,Med+3*Mad)) &  E[,y2[i]]<=Med+30*Mad  ))
#      }
#    }
   rm1<-which(V1/2<=35)
   rm2<-which(V2/2<=35)
   if(length(rm1)!=0){
     x3<-x2[-rm1]
     V1.3<-V1[-rm1]
   }else{
     x3<-x2
     V1.3<-V1
   }
   if(length(rm2)!=0){
     y3<-y2[-rm2]
     V2.3<-V2[-rm2]
   }else{
     y3<-y2
     V2.3<-V2
   }
x3M<-mean(x3)
y3M<-mean(y3)

x3S<-sd(x3)
y3S<-sd(y3)

  if(dim(E)[1]<900 & dim(E)[2]<900){
    dx3<-x3[-length(x3)]
    dV1.3<-0
    dy3<-y3[-length(y3)]
    dV2.3<-0
    for(k1 in 1:(length(x3)-1)){
      dV1.3[k1]<- (V1.3[k1+1]-V1.3[k1])/(x3[k1+1]-x3[k1]) 
    }
    for(k2 in 1:(length(y3)-1)){
      dV2.3[k2]<- (V2.3[k2+1]-V2.3[k2])/(y3[k2+1]-y3[k2]) 
    }
    dV2.3[length(V2.3)]<-dV2.3[length(dV2.3)]
    dV1.3[length(V1.3)]<-dV1.3[length(dV1.3)]
    nxd<-which(abs(dV1.3)<=1)
    nyd<-which(abs(dV2.3)<=1)
    V2.3<-V2.3[-nyd]
    y3<-y3[-nyd]
    V1.3<-V1.3[-nxd]
    x3<-x3[-nxd]
  }

     f<-function(p,x,V1){
      if(p[1]>0 & p[1]<max(x) & p[3]>0){
       r<-286.1
       if( (length(x)>10)){
         x0<-x
         y0<-V1/p[3]
         alpha<-atan2(y0,x0-p[1])
         Dis<- sqrt( (  x0 - p[2]*cos(alpha) - p[1])**2 + ( y0 - p[2]*sin(alpha) )**2 )
         rem<-which(Dis>20)
          if(length(rem)!=0){
             dis<-Dis[-rem]
           }
        return(sum( Dis  ))
       }else{return(1e300)}
      }else{return(1e300)}
     }

     f0<-function(p,x,V1){
       r<-286.1
       x0<-x
       y0<-V1/p[3]
       alpha<-atan2(y0,x0-p[1])
       Ms<-list(x= p[2]*cos(alpha) + p[1],y= p[2]*sin(alpha)*p[3])
       return(Ms)
     }


     f2<-function(p,x,V1){
      if(p[1]>0 & p[1]<max(x) & p[3]>0){
       r<-286.1
       neg<-which(  p[2]**2 < (x-p[1])**2 )
       if(length(neg)!=0){
          x<-x[-neg]
       }
       if(length(x)>10){
         y<-sqrt( p[2]**2 - (x-p[1])**2  )
         NAs<-which(is.na(y)==TRUE)
         if( length(NAs)!=0 & (length(x)-length(NAs))>10){
          return(sum( (y[-NAs]-V1[-NAs]/p[3])**2  ))
         }else{
          return(sum( (y-V1/p[3])**2  ))
         }
       }else{return(1e200)}
      }else{return(1e200)}
     }

   P0x<-matrix(nr=0,nc=3)
   P0x<-rbind(P0x,c(x3M,287,1.7))
   P0x<-rbind(P0x,c(x3M,290,1.7))
   P0x<-rbind(P0x,c(x3M,287,2.0))
   P0x<-rbind(P0x,c(x3M,287,1.5))
   P0x<-rbind(P0x,c(x3M,286,1.7))
   P0x<-rbind(P0x,c(x3M,288,1.7))
   P0x<-rbind(P0x,c(x3M,285,2))
   for(ip in -4:4){
     P0x<-rbind(P0x,c(x3M+i*x3S,287,1.7))
   }

   P0y<-matrix(nr=0,nc=3)
   P0y<-rbind(P0y,c(y3M,287,1.7))
   P0y<-rbind(P0y,c(y3M,290,1.7))
   P0y<-rbind(P0y,c(y3M,287,2.0))
   P0y<-rbind(P0y,c(y3M,287,1.5))
   P0y<-rbind(P0y,c(y3M,286,1.7))
   P0y<-rbind(P0y,c(y3M,288,1.7))
   P0y<-rbind(P0y,c(y3M,285,2))
   for(ip in -4:4){
     P0y<-rbind(P0y,c(y3M+i*y3S,287,1.7))
   }


   IsOk<-1
   eiX<-1
   eiY<-1
   while(IsOk!=0){
    if(class(try(nlm(f,p=P0x[eiX,],x=x3,V1=V1.3,hessian=TRUE,iterlim=8000,stepmax=100,steptol=1e-7)))!='try-error'){
      Xap<-nlm(f,p=P0x[eiX,],x=x3,V1=V1.3,hessian=TRUE,iterlim=8000,stepmax=100,steptol=1e-7)
      IsOk<-0      
    }else{
      if(eiX!=length(P0x[,1])){
        eiX<-eiX+1
      }else{
       Xap<-list(result=c(9999,9999,9999),min=0)
       IsOk<-0
      }
    }
   }
   IsOk<-1
   while(IsOk!=0){
    if(class(try(nlm(f,p=P0y[eiY,],x=y3,V1=V2.3,hessian=TRUE,iterlim=8000,stepmax=100,steptol=1e-7)))!='try-error'){
      Yap<-nlm(f,p=P0y[eiY,],x=y3,V1=V2.3,hessian=TRUE,iterlim=8000,stepmax=100,steptol=1e-7)
      IsOk<-0      
    }else{
      if(eiY!=length(P0y[,1])){
        eiY<-eiY+1
      }else{
        Yap<-list(result=c(9999,9999,9999),min=0)
        IsOk<-0
      }
    }
   }

###
#   aplot(x3,V1.3)
#   lines(f0(Xap$est,x3,V1.3),col="red")

#   aplot(y3,V2.3)
#   lines(f0(Yap$est,y3,V2.3),col="red")
###
    if(Xap$min!=0){
     if(class(try(sqrt(abs(diag(3*Xap$minimum/(length(x3) - 3)*solve(Xap$hessian))))))!='try-error'){
       ErrX<-sqrt(abs(diag(3*Xap$minimum/(length(x3) - 3)*solve(Xap$hessian))))
     }else{
       ErrX<-c(0,0,0)
     }
    }else{
     ErrX<-c(0,0,0)
    }
    if(Yap$min!=0){
     if(class(try(sqrt(abs(diag(3*Yap$minimum/(length(y3) - 3)*solve(Yap$hessian))))))!='try-error'){
       ErrY<-sqrt(abs(diag(3*Yap$minimum/(length(y3) - 3)*solve(Yap$hessian))))
     }else{
       ErrY<-c(0,0,0)
     }
    }else{
     ErrY<-c(0,0,0)
    }   
     Xap0<-c(Xap$est[1],ErrX[1]*1.96)  
     Yap0<-c(Yap$est[1],ErrY[1]*1.96)

     return(c(Xap0,Yap0))

    }
    x<-1:dim(E)[1]
    y<-1:dim(E)[2]
    cent<-XYap2(E,x,y,1)

    return(cent) # return the Xc Rx Yc Ry
  }

M<-matrix(nr=length(names),nc=9)
  Test<-read.table(paste0(work.area,"Super_F_table.txt"),header=TRUE)
for(i in 1:length(names)){
#for(i in 1:2){
  cat(paste0("1st  --- ",i," / ",length(names)," ...\n"))
   Img<-readFITS(paste0(path.data,as.character(names[i])))$imDat
###
  aux1aux<-which(paste0("Scs_",Test$FILE)==as.character(names[i]))
  aux2aux<-which(paste0("Sky_",Test$FILE)==as.character(names[i]))
  if(length(aux1aux)!=0){
    IND<-aux1aux
  }else{
    IND<-aux2aux
  }

###
   Xc<- round(Test$XcenC[IND])
   Yc<- round(Test$YcenC[IND])
   R<-40*2
   Z<-Img[(Xc-R):(Xc+R),(Yc-R):(Yc+R)]
###

#   png("PP.png")
#   image.plot(Z,zlim=c( median(Z) + 5*mad(Z)  ,max(Z)))
#   dev.off()

   indexs<-which(Z>median(Z,na.rm=T) + 5*mad(Z,na.rm=T),arr.ind=T)
  if( Xc!=-9999  ){
   Xc<-((Xc-R):(Xc+R))[round(median(indexs[,1],na.rm=T))]   
   Yc<-((Yc-R):(Yc+R))[round(median(indexs[,2],na.rm=T))]
   R<-40
   Z<-Img[(Xc-R):(Xc+R),(Yc-R):(Yc+R)]

   HotPix<-function(Z){
    Xaux<-1:dim(Z)[1]
    Yaux<-1:dim(Z)[2]
    XYaux<-expand.grid(Xaux,Yaux)
    Aux0<-which(  sqrt(  (XYaux[,1]-mean(Xaux))**2 + (XYaux[,2]-mean(Yaux))**2  )<=17  )
    Zaux00<-Z
    Zaux00[Aux0]<-NA
    aaaa<-which(  abs(Zaux00-mean(Zaux00,na.rm=T) )>= 10 *sd(Zaux00,na.rm=T)  ,arr.ind=T)
  while(length(aaaa)>0){
    Zaux<-Z
    Zaux[Aux0]<-NA
    aaaa<-which(  abs(Zaux-mean(Zaux,na.rm=T) )>= 10 *sd(Zaux,na.rm=T)  ,arr.ind=T)
#   if(length(aaaa)!=0){
#
#    image.plot(1:dim(Z)[1],1:dim(Z)[2],Zaux)
#
#    for(aaak00 in 1:length(aaaa[,1])){points(aaaa[aaak00,1],aaaa[aaak00,2])}   
#    aakk<-interp.surface(list(x=Xaux,y=Yaux,z=Z),aaaa)
   
    Zaux2<-Z
    for(aaak in 1:length(aaaa[,1])){
      Zaux2[aaaa[aaak,1],aaaa[aaak,2]]<-NA
    } 

   Rem.NAs<-function(M,OLD){
    Mc<-M
    UT<-1
    UT.aux<-0
    aux<-0
    if(OLD=="TRUE"){
      Mc[1:2,]<-0
      Mc[,1:2]<-0
      Mc[dim(Mc)[1]-(0:1),]<-0
      Mc[,dim(Mc)[2]-(0:1)]<-0
    }
    while(UT!=0){
      library(zoo)
      Na1<-length(which(is.na(Mc)==TRUE))
      Mc<-na.approx(Mc,na.rm = FALSE)
      Na2<-length(which(is.na(Mc)==TRUE)) 
      if(Na1==Na2){
        Mc<-t(na.approx(t(Mc),na.rm = FALSE))
      }
      UT<-length(which(is.na(Mc)==TRUE)) 
    }
    return(Mc)
   }
    Znew<-Rem.NAs(Zaux2,FALSE)
    Z<-Znew  #    aaaa<-which(  abs(Z-mean(Z,na.rm=T) )>= 10 *sd(Z,na.rm=T)  ,arr.ind=T)

    Zaux<-Z
    Zaux[Aux0]<-NA
    aaaa<-which(  abs(Zaux-mean(Zaux,na.rm=T) )>= 10 *sd(Zaux,na.rm=T)  ,arr.ind=T)
#   }
   }
   return(Z)
  }
  Z<-HotPix(Z)
  }


###
#png("KK.png");image.plot(Z);dev.off()
###
   Z<-Z/max(Z)
   Z<-Z-median(Z,na.rm=T)
   Z<-Z/sum(Z)
   X<-(Xc-R):(Xc+R)
   Y<-(Yc-R):(Yc+R)
   Zx<-0
   for(iq in 1:length(X)){
    Zx[iq]<-sum(Z[iq,])
   }
   Zy<-0
   for(iu in 1:length(Y)){
    Zy[iu]<-sum(Z[,iu])
   }

   nonX<-which(Zx<=0.015)
   nonY<-which(Zy<=0.015)

   if(length(nonX)!=0){
    X<-X[-nonX]
    Zx<-Zx[-nonX]
   }
   if(length(nonY)!=0){
    Y<-Y[-nonY]
    Zy<-Zy[-nonY]
   }
  
   Pn<-matrix(nr=0,nc=4)
   Pn<-rbind(Pn,c(mean(X),mean(Y),5,5))
   Pn<-rbind(Pn,c(mean(X)+sd(X),mean(Y),5,5))
   Pn<-rbind(Pn,c(mean(X)-sd(X),mean(Y),5,5))
   Pn<-rbind(Pn,c(mean(X),mean(Y)+sd(Y),5,5))
   Pn<-rbind(Pn,c(mean(X),mean(Y)-sd(Y),5,5))
   Pn<-rbind(Pn,c(mean(X)+sd(X),mean(Y)+sd(Y),5,5))
   Pn<-rbind(Pn,c(mean(X)-sd(X),mean(Y)-sd(Y),5,5))
   Pn<-rbind(Pn,c(mean(X),mean(Y),4,4))
   Pn<-rbind(Pn,c(mean(X),mean(Y),6,6))
   Pn<-rbind(Pn,c(mean(X),mean(Y),4,6))
   xx1<-which(Zx==max(Zx))
   yy1<-which(Zy==max(Zy))

  PNx<-matrix(nr=0,nc=6)
  PNx<-rbind(PNx,c(X[xx1],15,0.04,2.5,0.05,0.01))
  PNx<-rbind(PNx,c(X[xx1],10,0.04,1.7,0.03,0.05))
  PNx<-rbind(PNx,c(X[xx1],17,0.06,2.1,0.06,0.03))
  PNx<-rbind(PNx,c(X[xx1],13,0.04,2.5,0.05,0.03))
  PNx<-rbind(PNx,c(X[xx1],15,0.04,2.5,0.05,0.01))
  PNx<-rbind(PNx,c(mean(X),15,0.04,2.5,0.05,0.01))

  PNy<-matrix(nr=0,nc=6)
  PNy<-rbind(PNy,c(Y[yy1],15,0.04,2.5,0.05,0.01))
  PNy<-rbind(PNy,c(Y[yy1],10,0.04,1.7,0.03,0.05))
  PNy<-rbind(PNy,c(Y[yy1],17,0.06,2.1,0.06,0.03))
  PNy<-rbind(PNy,c(Y[yy1],13,0.04,2.5,0.05,0.03))
  PNy<-rbind(PNy,c(Y[yy1],15,0.04,2.5,0.05,0.01))
  PNy<-rbind(PNy,c(mean(Y),15,0.04,2.5,0.05,0.01))

  IsOkX<-1
  eiX<-1
  while(IsOkX!=0){
   if(class(try(resX<-nlm(Fun2,p=PNx[eiX,],X=X,Z=Zx,hessian=TRUE,
      iterlim=8000,steptol = 1e-8,gradtol = 1e-9,stepmax=1000)))!='try-error'){
    resX<-nlm(Fun2,p=PNx[eiX,],X=X,Z=Zx,hessian=TRUE,iterlim=8000,steptol = 1e-8,gradtol = 1e-9,stepmax=1000)
    IsOkX<-0
   }else{
     if(eiX!=length(PNx[,1])){
       eiX<-eiX+1
     }else{
      resX<-list(est=c(9999,9999,9999,9999),min=999999,hessian=matrix(0,nr=6,nc=6))
      IsOkX<-0
     }
   }
   }

  IsOkY<-1
  eiY<-1
  while(IsOkY!=0){
   if(class(try(resY<-nlm(Fun2,p=PNy[eiY,],X=Y,Z=Zy,hessian=TRUE,
      iterlim=8000,steptol = 1e-8,gradtol = 1e-9,stepmax=1000)))!='try-error'){
    resY<-nlm(Fun2,p=PNy[eiY,],X=Y,Z=Zy,hessian=TRUE,iterlim=8000,steptol = 1e-8,gradtol = 1e-9,stepmax=1000)
    IsOkY<-0
   }else{
     if(eiY!=length(PNy[,1])){
       eiY<-eiY+1
     }else{
      resY<-list(est=c(9999,9999,9999,9999),min=999999,hessian=matrix(0,nr=6,nc=6))
      IsOkY<-0
     }
   }
   }

#
system(paste0("mkdir ",work.area,"/Plots_kk"))
if(Xc!=-9999){
#png(paste0(work.area,"/Plots_kk/Fit_kk_XY",i,".png"))
#par(mfrow=c(1,2))
# aplot(X,Zx,type="l") 
#   kkfun<-function(p,X){
#     y<-p[3]*exp( -((X-p[1])**2)/(2*p[2])  ) + p[6]/2
#     y2<-p[5]*exp( -((X-p[1])**2)/(2*p[4])  ) + p[6]/2
#     S<-y+y2
#   return(S)}
# lines(X,kkfun(resX$est,X),col="red")
#
# aplot(Y,Zy,type="l") 
# lines(Y,kkfun(resY$est,Y),col="blue")
#dev.off()
}
#

   if(class(try(sqrt(abs(diag(6*resX$minimum/(length(X) - 6)*solve(resX$hessian))))))!='try-error'){
      ErrX<-sqrt(abs(diag(6*resX$minimum/(length(X) - 6)*solve(resX$hessian))))
   }else{
      ErrX<-c(0,0,0,0)
   }

   if(class(try(sqrt(abs(diag(6*resY$minimum/(length(Y) - 6)*solve(resY$hessian))))))!='try-error'){
      ErrY<-sqrt(abs(diag(6*resY$minimum/(length(Y) - 6)*solve(resY$hessian))))
   }else{
      ErrY<-c(0,0,0,0)
   }

   Cen<-CENT2(Img)

   M[i,2]<- resX$est[1]
   M[i,3]<- ErrX[1]*1.96
   M[i,4]<- resY$est[1]
   M[i,5]<- ErrY[1]*1.96

   M[i,6]<-Cen[1]
   M[i,7]<-Cen[2]
   M[i,8]<-Cen[3]
   M[i,9]<-Cen[4]
  cat(paste0("2d  --- ",i," / ",length(names)," ...\n"))
  cat(" ....\n ")
}



M[,1]<-names
colnames(M)<-c("Name","X","Err_X","Y","Err_Y","Xc","Err_Xc","Yc","Err_Yc")

write.table(M,paste0(work.area,"/Center_ScsSky_new3.txt"),row.names=FALSE,col.names=TRUE,quote=FALSE)


A<-read.table("Super_F_table.txt",header=TRUE)


##########

KKsc<-gsub("Scs_","",M[,1])
KKskc<-gsub("Sky_","",KKsc)

KKflag<-0
KKindex<-0
KKindex2<-0
namesSK<-KKskc
for(i in 1:length(namesSK)){
  auxKK<-which(namesSK[i]==A$FILE)
  KKindex[i]<-i
  KKindex2[i]<-auxKK[1]
  if(as.character(A$DPR.TYPE[auxKK[1]])=="SKY"){
   KKflag[i]<-2
  }else if(as.character(A$DPR.TYPE[auxKK[1]])=="OBJECT"){
   KKflag[i]<-1
  }else{
   KKflag[i]<-0
  }
}

Sky_im<-which(KKflag==2)

NSky<-cbind(A[KKindex2[Sky_im],],M[KKindex[Sky_im],])
NScs<-M[KKindex[-Sky_im],]
##########



#N<-matrix(nr=0,nc=dim(A)[2]+dim(M)[2])
#ValSky<-0
#ValSkyM<-0
#for(i in 1:length(names)){
# valSky<-which(paste0("Sky_",A$FILE)==M[i,1])
# ValSky[i]<-valSky[1]
# if(length(valSky)!=0){
#   ValSkyM[i]<-i
# }else{
#   ValSkyM[i]<-0
# }
#}

#ValScs<-0
#ValScsM<-0
#for(i in 1:length(names)){
# valScs<-which(paste0("Scs_",A$FILE)==M[i,1])
# if(length(valScs)!=0){
#   ValScsM[i]<-i
# }else{
#   ValScsM[i]<-0
# }
#}

#NASky<-which(is.na(ValSky)==TRUE)
#if(length(NASky)!=0){
#  ValSky<-ValSky[-NASky]
#  ValSkyM<-ValSkyM[-NASky]
#}

#NSky<-cbind(A[ValSky,],M[ValSkyM,])
#NScs<-M[ValScs,]

write.table(NSky,paste0(work.area,"/Match_Sky_center_new2.txt"),row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(NScs,paste0(work.area,"/Center_Scs_new2.txt"),row.names=FALSE,col.names=TRUE,quote=FALSE)

KK<-read.table(paste0(work.area,"/Match_Sky_center_new2.txt"),header=TRUE)

png(paste0(work.area,"/Centering_trans_corr.png"))

XL<-range(c(KK$X-KK$X[1],KK$Xc-KK$Xc[1]))
YL<-range(c(KK$Y-KK$Y[1],KK$Yc-KK$Yc[1]))
aplot(KK$X-KK$X[1],KK$Y-KK$Y[1],pch=20,col="red",xlab="X - X[1]",ylab="Y - Y[1]",xlim=XL,ylim=YL)
for(ip in 1:length(KK$X)){
  lines( c(KK$X[ip]-KK$X[1],KK$X[ip]-KK$X[1]) , c( KK$Y[ip]-KK$Y[1] - KK$Err_Y[ip],KK$Y[ip]-KK$Y[1] + KK$Err_Y[ip]  )  )
  lines( c(KK$X[ip]-KK$X[1]-KK$Err_X[ip],KK$X[ip]-KK$X[1]+KK$Err_X[ip]) , c( KK$Y[ip]-KK$Y[1],KK$Y[ip]-KK$Y[1]  )  )
  lines( c(KK$Xc[ip]-KK$Xc[1],KK$Xc[ip]-KK$Xc[1]) , c( KK$Yc[ip]-KK$Yc[1] - KK$Err_Yc[ip],KK$Yc[ip]-KK$Yc[1] + KK$Err_Yc[ip]  )  )
  lines( c(KK$Xc[ip]-KK$Xc[1]-KK$Err_Xc[ip],KK$Xc[ip]-KK$Xc[1]+KK$Err_Xc[ip]) , c( KK$Yc[ip]-KK$Yc[1],KK$Yc[ip]-KK$Yc[1]  )  )
}
points(KK$X-KK$X[1],KK$Y-KK$Y[1],pch=20,col="red")
points(KK$Xc-KK$Xc[1],KK$Yc-KK$Yc[1],pch=20,col="blue")
dev.off()


png(paste0(work.area,"/Centering_trans_corr_range.png"))

XL<-  c(median(c(KK$X-KK$X[1],KK$Xc-KK$Xc[1])) - 5*mad(c(KK$X-KK$X[1],KK$Xc-KK$Xc[1])),median(c(KK$X-KK$X[1],KK$Xc-KK$Xc[1])) + 5*mad(c(KK$X-KK$X[1],KK$Xc-KK$Xc[1])))
YL<-  c(median(c(KK$Y-KK$Y[1],KK$Yc-KK$Yc[1])) - 5*mad(c(KK$Y-KK$Y[1],KK$Yc-KK$Yc[1])),median(c(KK$Y-KK$Y[1],KK$Yc-KK$Yc[1])) + 5*mad(c(KK$Y-KK$Y[1],KK$Yc-KK$Yc[1])))
aplot(KK$X-KK$X[1],KK$Y-KK$Y[1],pch=20,col="red",xlab="X - X[1]",ylab="Y - Y[1]",xlim=XL,ylim=YL)
for(ip in 1:length(KK$X)){
  lines( c(KK$X[ip]-KK$X[1],KK$X[ip]-KK$X[1]) , c( KK$Y[ip]-KK$Y[1] - KK$Err_Y[ip],KK$Y[ip]-KK$Y[1] + KK$Err_Y[ip]  )  )
  lines( c(KK$X[ip]-KK$X[1]-KK$Err_X[ip],KK$X[ip]-KK$X[1]+KK$Err_X[ip]) , c( KK$Y[ip]-KK$Y[1],KK$Y[ip]-KK$Y[1]  )  )
  lines( c(KK$Xc[ip]-KK$Xc[1],KK$Xc[ip]-KK$Xc[1]) , c( KK$Yc[ip]-KK$Yc[1] - KK$Err_Yc[ip],KK$Yc[ip]-KK$Yc[1] + KK$Err_Yc[ip]  )  )
  lines( c(KK$Xc[ip]-KK$Xc[1]-KK$Err_Xc[ip],KK$Xc[ip]-KK$Xc[1]+KK$Err_Xc[ip]) , c( KK$Yc[ip]-KK$Yc[1],KK$Yc[ip]-KK$Yc[1]  )  )
}
points(KK$X-KK$X[1],KK$Y-KK$Y[1],pch=20,col="red")
points(KK$Xc-KK$Xc[1],KK$Yc-KK$Yc[1],pch=20,col="blue")
dev.off()





KKsk<-read.table(paste0(work.area,"/Match_Sky_center_new2.txt"),header=TRUE)

KKsc<-read.table(paste0(work.area,"/Center_Scs_new2.txt"),header=TRUE)

png(paste0(work.area,"/Big_circle.png"))

aplot(KKsc$Xc,KKsc$Yc,pch=20,col="blue")
points(KKsk$Xc,KKsk$Yc,pch=20,col="red")

dev.off()






