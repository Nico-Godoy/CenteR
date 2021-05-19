
REDUCTION<-function(work.area,path.data,NODES,OLD,...){
   setwd(path.data)
   system(paste0("dfits *.fits | fitsort DPR.CATG DPR.TYPE MJD-OBS EXPTIME NAXIS1 NAXIS2 \"HIERARCH ESO TPL NEXP\"> ",work.area,"/Info.list"))
   setwd(work.area)
   system(paste0("sed -i \'s|HIERARCH ESO TPL NEXP|NEXP|g\' ",work.area,"/Info.list"))
   R<-read.table(paste0(work.area,"/Info.list"),header=TRUE)

   sc<-which(R$DPR.CATG=="SCIENCE") 

   timesSc<-R$EXPTIME[sc]
   TimeSc<-sort(unique(R$EXPTIME[sc]),decreasing=TRUE)   
   if(length(TimeSc)>2){
    auxTime<-0
    for(i in 1:length(TimeSc)){
     auxTime[i]<-length(which(timesSc==TimeSc[i]))
    }
    time1<-which(auxTime>10)
#    time2<-which(auxTime==6 | auxTime==3 | auxTime==2)
     time2<-which(auxTime==6 | auxTime==2 | auxTime==4)
    timeOut<-which(auxTime==1)
    TimeSc<-c(TimeSc[time1],TimeSc[time2])
   }

   sc1<-which(R$EXPTIME==TimeSc[1] & R$DPR.CATG=="SCIENCE")
   sc2<-which(R$EXPTIME==TimeSc[2] & R$DPR.CATG=="SCIENCE")

   if(length(unique(R$NAXIS1[sc2]))!=1){
      dimsX<-unique(R$NAXIS1[sc2])
      lenX<-0
      for(lenii in 1:length(dimsX)){
        lenX[lenii]<-length(which(R$NAXIS1[sc2]==dimsX[lenii]))
      }
      rightLen<-which(R$NAXIS1[sc2]==dimsX[which(lenX==max(lenX))])
      sc2<-sc2[rightLen]
   }

   fl<-which(R$DPR.CATG=="CALIB" & R$DPR.TYPE=="FLAT,SKY")
   TimeFl<-unique(R$EXPTIME[fl])
   if(length(TimeFl)>1){
     TimeFl<-TimeFl[round(length(TimeFl)/2)]
   }
   fl1<-which(R$EXPTIME==TimeFl[1] & R$DPR.CATG=="CALIB" & R$DPR.TYPE=="FLAT,SKY")

  scC<-which(R$DPR.CATG=="CALIB" & R$DPR.TYPE=="DARK" & R$EXPTIME==TimeSc[1] & R$NAXIS1[sc1[1]]==R$NAXIS1 & R$NAXIS2[sc1[1]]==R$NAXIS2)  # Calib
  scS<-which(R$DPR.CATG=="CALIB" & R$DPR.TYPE=="DARK" & R$EXPTIME==TimeSc[2] & R$NAXIS1[sc2[1]]==R$NAXIS1 & R$NAXIS2[sc2[1]]==R$NAXIS2) # Sc
  flt<-which(R$DPR.CATG=="CALIB" & R$DPR.TYPE=="DARK" & R$EXPTIME==TimeFl[1] & R$NAXIS1[fl1[1]]==R$NAXIS1 & R$NAXIS2[fl1[1]]==R$NAXIS2) # Sc

  NamesscC<-R$FILE[scC]
  NamesscS<-R$FILE[scS]
  Namesflt<-R$FILE[flt]

  DimsC<-c(R$NAXIS1[sc1[1]],R$NAXIS2[sc1[1]])
  DimsS<-c(R$NAXIS1[sc2[1]],R$NAXIS2[sc2[1]])
  DimsF<-c(R$NAXIS1[fl1[1]],R$NAXIS2[fl1[1]])
 
  Master.Dark<-function(list,Dims){
    M<-array(0,dim=c(Dims,0))
    N<-matrix(0,nr=Dims[1],nc=Dims[2])
    for(i in 1:length(list)){
      K<-readFITS(paste0(path.data,"/",as.character(list[i])))$imDat
      M<-abind(M,K)
    }
    for(i in 1:dim(M)[1]){
     for(j in 1:dim(M)[2]){
      N[i,j]<-mean(M[i,j,],na.rm=TRUE)
     }}
    return(N)
  }

  MD.C<-Master.Dark(NamesscC,DimsC)
if(is.na(DimsS)!=TRUE){  MD.S<-Master.Dark(NamesscS,DimsS)}
  MD.f<-Master.Dark(Namesflt,DimsF)

####
  BadPixels<-function(A,Sigma,OLD){
     aw<-which(abs(A-median(A))>=Sigma*mad(A))
#     aw<-which(A>100)
     B<-A*0
     if(length(aw)>0){
       B[aw]<-10000
     }else{
       B[1,1]<-10000
     }
    if(OLD=="TRUE"){
     B[1:2,]<-0
     B[,1:2]<-0
     B[dim(B)[1]-(0:1),]<-0
     B[,dim(B)[2]-(0:1)]<-0
    }
    return(B)
  }

  BP_MD.C<-BadPixels(MD.C,5,OLD)
if(is.na(DimsS)!=TRUE){BP_MD.S<-BadPixels(MD.S,5,OLD)}
####

  system(paste0("mkdir ",work.area,"/Data"))
  writeFITSim(MD.C,paste0(work.area,"/Data/MD_",TimeSc[1],"_.fits"))
if(is.na(DimsS)!=TRUE){  writeFITSim(MD.S,paste0(work.area,"/Data/MD_",TimeSc[2],"_.fits"))}
  writeFITSim(MD.f,paste0(work.area,"/Data/MD_",TimeFl,"_.fits"))

####
  writeFITSim(BP_MD.C,paste0(work.area,"/Data/MD_",TimeSc[1],"_BP.fits"))
if(is.na(DimsS)!=TRUE){  writeFITSim(BP_MD.S,paste0(work.area,"/Data/MD_",TimeSc[2],"_BP.fits"))}
####

 Master.Flat<-function(fl,MD.f,Dims){
  FL<-array(dim=c(Dims,0))
  for(i in 1:length(fl)){
    K<-readFITS(paste0(path.data,"/",as.character(R$FILE[fl[i]])))$imDat
    dimZ<-dim(K)[3]
    if(is.na(dimZ)==FALSE){
      for(j in 1:dimZ){ FL<-abind(FL,K[,,j]-MD.f)}
    }else{  FL<-abind(FL,K-MD.f) }
  }
  MFL<-matrix(nr=Dims[1],nc=Dims[2])
  for(i in 1:dim(FL)[1]){
   for(j in 1:dim(FL)[2]){
    MFL[i,j]<-mean(FL[i,j,])
   }
  }
  MFL<-MFL/median(MFL)
  return(MFL)
 }

   MF.f<-Master.Flat(fl1,MD.f,DimsF)

   writeFITSim(MF.f,paste0(work.area,"/Data/MF.fits"))

   C1<-((DimsF[1]-DimsC[1])/2+1):(DimsC[1] + ((DimsF[1]-DimsC[1])/2))
   C2<-((DimsF[2]-DimsC[2])/2+1):(DimsC[2] + ((DimsF[2]-DimsC[2])/2))

####
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
####
   cl <- makeCluster(NODES)
   registerDoParallel(cl)  
   KK<-foreach(i=1:length(sc1),.packages=c('foreach','abind','FITSio','fields','akima','h5','astro'),.export=c("R","C1","C2","path.data","MD.C","MF.f","BP_MD.C","work.area"),.inorder=TRUE)%dopar%{
#   for(i in 1:length(sc1)){
     E<-readFITS(paste0(path.data,"/",as.character(R$FILE[sc1[i]])))$imDat
####
     NasC<-which(BP_MD.C>0)
####
     for(j in 1:dim(E)[3]){
      E[,,j]<-(E[,,j]-MD.C)/MF.f[C1,C2]
      E0<-E[,,j]
      E0[NasC]<-NA
      E0<-Rem.NAs(E0,OLD)

      if(OLD=="TRUE"){
       NasCK<-which(E0<= -1000)
       if(length(NasCK)!=0){
          E0[NasCK]<-NA
          E0<-Rem.NAs(E0,OLD)
       }
       NasCK0<-which(E0>=25000)
       if(length(NasCK)!=0){
          E0[NasCK0]<-NA
          E0<-Rem.NAs(E0,OLD)
       }
       E0.BP<-BadPixels(E0,5,OLD)
       NasCK1<-which(BP_MD.C>0)
       if(length(NasCK1)!=0){
          E0[NasCK1]<-NA
          E0<-Rem.NAs(E0,OLD)
       }
      }
      E[,,j]<-E0
     }
     infg<-which(E==Inf | is.nan(E)==TRUE | is.na(E)==TRUE)
     if(length(infg)!=0){E[infg]<-0}
     writeFITSim(E,paste0(work.area,"/Data/",gsub(".fits","",R$FILE[sc1[i]]),"_DF.fits"))
   }
   stopCluster(cl)

if(is.na(DimsS)!=TRUE){
   cl <- makeCluster(NODES)
   registerDoParallel(cl)  
   S1<-((DimsF[1]-DimsS[1])/2+1):(DimsS[1] + ((DimsF[1]-DimsS[1])/2))
   S2<-((DimsF[2]-DimsS[2])/2+1):(DimsS[2] + ((DimsF[2]-DimsS[2])/2))

   KK<-foreach(i=1:length(sc2),.packages=c('foreach','abind','FITSio','fields','akima','h5','astro'),.export=c("R","S1","S2","path.data","MD.S","MF.f","BP_MD.S","work.area"),.inorder=TRUE)%dopar%{
#   for(i in 1:length(sc2)){
     E<-readFITS(paste0(path.data,"/",as.character(R$FILE[sc2[i]])))$imDat
####
     NasS<-which(BP_MD.S>0)
####
     for(j in 1:dim(E)[3]){
cat(paste0(c(dim(MD.S),"  ---  ",dim(E),"\n")))
      E[,,j]<-(E[,,j]-MD.S)/MF.f[S1,S2]
      E0<-E[,,j]
      E0[NasS]<-NA
      E0<-Rem.NAs(E0,OLD)

      if(OLD=="TRUE"){
       NasCK<-which(E0<= -1000)
       if(length(NasCK)!=0){
          E0[NasCK]<-NA
          E0<-Rem.NAs(E0,OLD)
       }
       NasCK0<-which(E0>=25000)
       if(length(NasCK0)!=0){
          E0[NasCK0]<-NA
          E0<-Rem.NAs(E0,OLD)
       }
       E0.BP<-BadPixels(E0,5,OLD)
       NasCK1<-which(BP_MD.S>0)
       if(length(NasCK1)!=0){
          E0[NasCK1]<-NA
          E0<-Rem.NAs(E0,OLD)
       }
      }
      E[,,j]<-E0
     }
     writeFITSim(E,paste0(work.area,"/Data/",gsub(".fits","",R$FILE[sc2[i]]),"_DF.fits"))
   }
   stopCluster(cl)
}
  write.table(R,paste0(work.area,"/R_val.txt"),col.names=TRUE,row.names=FALSE,quote=FALSE)
}

######
######
if(Data_reduction=="TRUE"){
 cat("Data reduction in process ...\n")
 REDUCTION(work.area2,path.data2,NODES,OLD)
}
######
######


#-------------------------------------------------------------------------------

Master<-function(path.data,work.area,ra,dec){

# path.data<-path.data1
# work.area<-work.area1
# ra<-ra1
# dec<-dec1

  cf<-180/pi
  RA<-hms2rad(ra)*cf
  DEC<-dms2rad(dec)*cf

##
  setwd(path.data)
  names<-system("ls *.fits",intern=TRUE)
  setwd(work.area)

  Max<-matrix(nr=length(names),nc=1) 
  for(i in 1:length(Max[,1])){
   Max[i,]<-paste0("dfits ",path.data,"/",names[i]," | fitsort RA DEC UTC \"HIERARCH ESO TEL AMBI FWHM START\" \"HIERARCH ESO TEL AMBI FWHM END\" \"HIERARCH ESO TEL IA FWHM\"  \"HIERARCH ESO TEL AMBI TAU0\" \"HIERARCH ESO TEL AMBI WINDSP\" \"HIERARCH ESO TEL AIRM START\" \"HIERARCH ESO TEL AIRM END\" \"HIERARCH ESO TEL FOCU VALUE\" NAXIS3  \"HIERARCH ESO TEL TARG DELTA\" \"HIERARCH ESO DPR TYPE\" \"HIERARCH ESO SEQ NEXPO\" > ",work.area,"/val")
  }
  Max2<-matrix(nr=length(Max[,1]),nc=14)
# 4
  Max2b<-matrix(nr=length(Max[,1]),nc=1)
#  Counts<-0
  for(i in 1:length(Max[,1])){
    system(Max[i,])
    kk<-read.table(paste0(work.area,"val"),sep="\t",header=TRUE)
    if(is.numeric(kk[1,2])=='FALSE'){   
       kk[1,2:12]<-999999999
    }
    Max2[i,2]<-kk[1,2]   # RA
    Max2[i,3]<-kk[1,3]   # DEC
    Max2[i,4]<-kk[1,4]   # Date/UTC
    Max2[i,5]<-kk[1,5]   # FWHM_start
    Max2[i,6]<-kk[1,6]   # FWHM_end
    Max2[i,7]<-kk[1,7]   # FWHM_IA
    Max2[i,8]<-kk[1,8]   # Coher_time
    Max2[i,9]<-kk[1,9]   # Wind_speed
    Max2[i,10]<-kk[1,10] # Airmass_start
    Max2[i,11]<-kk[1,11] # Airmass_end
    Max2[i,12]<-kk[1,12] # Focu_value
    Max2[i,13]<-kk[1,13] # z-Dim
    Max2[i,14]<-kk[1,16] # NEXPO
    Max2b[i,1]<-as.character(kk[1,14]) # DPOS
  }
  system(paste0("rm ",work.area,"/val"))

  Dist<-0

  for(i in 1:length(Max2[,1])){
   if(Max2[i,2]!=999999999){
     dist<-angSep(rad2hms(RA/cf,places=6),rad2dms(DEC/cf,places=6),rad2hms(Max2[i,2]/cf,places=6),rad2dms(Max2[i,3]/cf,places=6))
     Dist[i]<-dist*3600
   }else{
     Dist[i]<-999999999
   }
  }
  objt<-which(Dist<5)
  objt0<-which(Dist==999999999)  
  INtime<-sort.int(Max2[,4],index.return=TRUE)$ix

  Max2<-Max2[INtime,]
  Max2b<-Max2b[INtime,1]
  
  M1<-cbind(matrix(names[INtime],nc=1),Max2[,2:14],Dist[INtime])
 
  Header<-c(" Name RA DEC Time FWHM_start FWHM_end FWHM_IA Coher_time Wind_speed Airmass_start Airmass_end Focu_value DimZ NEXPO Dist")
  header<-matrix(Header,nr=1,nc=1)

  write.table(header,paste0(work.area,"/Master_inf.txt"),col.names=FALSE,row.names=FALSE,quote=FALSE,append=FALSE)
  write.table(M1    ,paste0(work.area,"/Master_inf.txt"),col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE)

}


######
######
if(Master_File=="TRUE"){
 cat("Creating the first Master file ...\n")
 Master(path.data2,work.area2,ra2,dec2)
}
######
######

#-------------------------------------------------------------------------------


CENTER<-function(){

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
       Med<-median(E[-NoU])
       Mad<-mad(E[-NoU])
       E[NoU]<- -99999
    }else{
       Med<-median(E)
       Mad<-mad(E)
    }
    for(i in 1:max(c(length(x2),length(y2)))){
      if(i<=length(x2)){
        V1[i]<-length(which( E[x2[i],]>=max(c(0,Med+3*Mad)) &  E[x2[i],]<=Med+30*Mad  ))
      }
      if(i<=length(y2)){
        V2[i]<-length(which( E[,y2[i]]>=max(c(0,Med+3*Mad)) &  E[,y2[i]]<=Med+30*Mad  ))
      }
    }
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
     f<-function(p,x,V1){
      if(p[1]>0 & p[1]<max(x)){
       r<-286.1
       if( (length(x)>10)){
         x0<-x
         y0<-V1/2
         alpha<-atan2(y0,x0-p[1])
         Dis<- sqrt( (  x0 - p[2]*cos(alpha) - p[1])**2 + ( y0 - p[2]*sin(alpha) )**2 )
         rem<-which(Dis>20)
          if(length(rem)!=0){
             dis<-Dis[-rem]
           }
        return(sum( Dis  ))
       }else{return(1e30)}
      }else{return(1e30)}
     }

     f2<-function(p,x,V1){
      if(p[1]>0 & p[1]<max(x)){
       r<-286.1
       y<-sqrt( p[2]**2 - (x-p[1])**2  )
       NAs<-which(is.na(y)==TRUE)
       if(length(NAs)!=0 & (length(x)-length(NAs))>10){
        return(sum( (y[-NAs]-V1[-NAs]/1.8)**2  ))
       }else{return(1e30)}
      }else{return(1e30)}
     }
    # added @ 25-03-2018
     NAs13<-which(is.na(x3)==TRUE | is.na(V1.3)==TRUE)
     if(length(NAs13)!=0){
       x3<-x3[-NAs13]
       V1.3<-V1.3[-NAs13]
     }
     NAs12<-which(is.na(y3)==TRUE | is.na(V2.3)==TRUE)
     if(length(NAs12)!=0){
       y3<-y3[-NAs12]
       V2.3<-V2.3[-NAs12]
     }
     if(length(x3)!=0){
        Xap<-nlm(f2,p=c(mean(x3),287),x=x3,V1=V1.3)$est
     }else{
        Xap<-c(-9999,-9999)
     }
     if(length(y3)!=0){
        Yap<-nlm(f2,p=c(mean(y3),287),x=y3,V1=V2.3)$est
     }else{
        Yap<-c(-9999,-9999)
     }
     return(c(Xap,Yap))
    }
    x<-1:dim(E)[1]
    y<-1:dim(E)[2]
    cent<-XYap2(E,x,y,2)
    return(cent) # return the Xc Rx Yc Ry
  }
###
  R<-read.table(paste0(work.area,"/R_val.txt"),header=TRUE)
  Aa<-read.table(paste0(work.area,"/Master_inf.txt"),header=TRUE)
  aas<-match(R$FILE,Aa$Name)

  A<-cbind(R,Aa[aas,-1])

   Ident<-function(Lista){
    Time<-range(Lista$Time)
    ok<-1
    while(ok==1){
     sc<-which(Lista$DPR.CATG=="SCIENCE") 
     TimeSc<-sort(unique(R$EXPTIME[sc]),decreasing=TRUE)
     sc10<-which(Lista$EXPTIME==TimeSc[1] & Lista$DPR.CATG=="SCIENCE" & Lista$Time>=Time[1] & Lista$Time<=Time[2] & is.na(Lista$DimZ)==FALSE)
     sc20<-which(Lista$EXPTIME==TimeSc[2] & Lista$DPR.CATG=="SCIENCE" & Lista$Time>=Time[1] & Lista$Time<=Time[2] & is.na(Lista$DimZ)==FALSE)
     if(length(sc10)>length(sc20)){
       sc1<-sc20
       sc2<-sc10
     }else{
       sc1<-sc10
       sc2<-sc20
     }
     List<-Lista[sc2,]
     poi<-matrix(nr=length(List$Time),nc=3)
     if( length(which(List$DPR.TYPE=='SKY'))!=0){
         sky<-which(List$DPR.TYPE=='SKY')
         scs<-which(List$DPR.TYPE=='OBJECT')
     }else{
         scs<-which(List$NEXPO>2)
         sky<-which(List$NEXPO<=2)
     }
     if(length(sky)<=4){ ## Added becuae a fail in the classification in new data
        KM0<-kmeans(List$Time,centers=2)$cluster

        Position<-1:length(List$Dist)

        KMa<-kmeans(List$Dist[which(KM0==1)],centers=2)$cluster
        C1a<-which(KMa==1)
        C2a<-which(KMa==2)
        if(length(C1a)>length(C2a)){
          scsa<-Position[which(KM0==1)][C1a]
          skya<-Position[which(KM0==1)][C2a]
        }else{
          scsa<-Position[which(KM0==1)][C2a]
          skya<-Position[which(KM0==1)][C1a]
        }

        KMb<-kmeans(List$Dist[which(KM0==2)],centers=2)$cluster
        C1b<-which(KMb==1)
        C2b<-which(KMb==2)
        if(length(C1b)>length(C2b)){
          scsb<-Position[which(KM0==2)][C1b]
          skyb<-Position[which(KM0==2)][C2b]
        }else{
          scsb<-Position[which(KM0==2)][C2b]
          skyb<-Position[which(KM0==2)][C1b]
        }
        scs<-c(scsa,scsb)
        sky<-c(skya,skyb)

      }
     ok<-2
    png(paste0(work.area,"/Classification_v2.png"),width = 680, height = 680)
       plot(List$Time,List$Dist,xlab="Time",ylab="Dist")
       points(List$Time[scs],List$Dist[scs],pch=20,col="blue")
       points(List$Time[sky],List$Dist[sky],pch=19,col="red")
    dev.off()
     }
     v0 <- dim(Lista)[2]
     LISTA<-cbind(Lista,matrix(0,nc=1,nr=dim(Lista)[1])) # pensar en una forma inteligente de ver esto
     names(LISTA)[v0+1]<-"ObjType"
     v1n<-match(List$FILE[scs],Lista$FILE)
     v2n<-match(List$FILE[sky],Lista$FILE)
     LISTA$ObjType[v1n]<-"Science"
     LISTA$ObjType[v2n]<-"Sky"
     Nans<-which(LISTA$ObjType=="0")
     if(length(Nans)!=0){
        LISTA$ObjType[Nans]<-"Nan"}
     return(LISTA)
   }

   List<-Ident(A)
   write.table(List,paste0(work.area,"/Master_info_v2.txt"),col.names=TRUE,row.names=FALSE,quote=FALSE)
   A<-read.table(paste0(work.area,"/Master_info_v2.txt"),header=TRUE)

  obj<-which(A$ObjType!="Nan")

  M<-matrix(nr=0,nc=dim(A)[2]+5) 
  for(i in 1:length(obj)){
   cat("----------------------------------\n")
   cat(paste0("---- ",i," / ",length(obj),": Fisrt part\n"))
   Im<-readFITS(paste0(work.area,"/Data/",   gsub(".fits","_DF.fits",as.character(A$FILE[obj[i]]))))$imDat
   if(A$ObjType[obj[i]]=="Science"){
     Im2<-Im[,,1]
     for(j2 in 1:dim(Im)[1]){
      for(k2 in 1:dim(Im)[2]){
       Im2[j2,k2]<-median(Im[j2,k2,])
      }
     }
     cat("-- prev...\n")
     cent<-CENT2(Im2)
     cat(paste0("---- ",i," / ",dim(A)[1],": Second part\n"))
        cat(paste0(" ---- center: ",cent[1]," ",cent[3],"\n"))
     M<-rbind(M,cbind(A[obj[i],],matrix(c( paste0(gsub(".fits","",as.character(A$FILE[obj[i]])),"_scs"),"Scs",1,cent[1],cent[3]),nr=1)))
     writeFITSim(Im2,paste0(work.area,"/Data/Scs_",as.character(A$FILE[obj[i]])))

   }else if(A$ObjType[obj[i]]=="Sky"){
     Im2<-Im[,,1]
     for(j2 in 1:dim(Im)[1]){
      for(k2 in 1:dim(Im)[2]){
        Im2[j2,k2]<-median(Im[j2,k2,])
      }
     }
     cat("-- prev...\n")
     cent<-CENT2(Im2)
     cat(paste0("---- ",i," / ",dim(A)[1],": Second part\n"))
        cat(paste0(" ---- center: ",cent[1]," ",cent[3],"\n"))
     M<-rbind(M,cbind(A[obj[i],],matrix(c( paste0(gsub(".fits","",as.character(A$FILE[obj[i]])),"_sky"),"Sky",1,cent[1],cent[3]),nr=1)))
     writeFITSim(Im2,paste0(work.area,"/Data/Sky_",as.character(A$FILE[obj[i]])))

   }else{
     M<-rbind(M,cbind(A[i,],c("Nam","noSc",0,0,0)))
   }
    cat(paste0("---- ",i," / ",dim(A)[1],": Third part\n"))
 }
 names(M)[dim(A)[2]+1]<-"NewName"
 names(M)[dim(A)[2]+2]<-"Flag"
 names(M)[dim(A)[2]+3]<-"numImage"
 names(M)[dim(A)[2]+4]<-"XcenC"
 names(M)[dim(A)[2]+5]<-"YcenC"
  cat("---- Writing the file ...\n")

  noPe<-which(M$XcenC<= -1)
  if(length(noPe)!=0){
    M<-M[-noPe,]
  }
 write.table(M,"Super_F_table.txt",col.names=TRUE,row.names=FALSE,quote=FALSE)

}

######
######
if(Center_File=="TRUE"){
 cat(" --- Running the Center file ...\n")
 CENTER()
}
######
######


write.table("ok1","ok1.txt",col.names=FALSE,row.names=FALSE,quote=FALSE)

###
SKYsub<-function(...){

RT<-read.table("Super_F_table.txt",header=TRUE)

O1<-which(RT$Flag=="Scs" | RT$Flag=="Science")
O2<-which(RT$Flag=="Sky")

Xc<-mean(RT$XcenC)
Yc<-mean(RT$YcenC)

RT<-read.table("Super_F_table.txt",header=TRUE)

SC1<-which(RT$Flag=="Scs" | RT$Flag=="Science")
SK1<-which(RT$Flag=="Sky")

for(i in 1:length(SC1)){
  cat(paste0(i,"/",length(SC1)," \n"))
  IMG<-readFITS(paste0(work.area,"/Data/",gsub(".fits","_DF.fits",as.character(RT$FILE[SC1[i]]))))$imDat
  tSk<-which(  abs(RT$Time[SK1]-RT$Time[SC1[i]])==min(abs(RT$Time[SK1]-RT$Time[SC1[i]])))
  SKY<-readFITS(paste0(work.area,"/Data/Sky_",as.character(RT$FILE[SK1[tSk[1]]])))$imDat
  for(j in 1:dim(IMG)[3]){
    IMG[,,j]<-IMG[,,j]-SKY
  }
  writeFITSim(IMG,paste0(work.area,"/Data/Scs_", gsub(".fits","_DF.fits",as.character(RT$FILE[SC1[i]])) ))
  rm(IMG)
}

write.table("ok2","ok1.txt",col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE)
}


######
######
if(Sky_Subs=="TRUE"){
 cat(" --- Running the Sky substraction ...\n")
 SKYsub()
}
###### 
######

library(IM)
  Test<-read.table("Super_F_table.txt",header=TRUE)

  ZNPcNC<-function(Ncomp,Name,Images){
   Xc<- median(Test$XcenC)
   Yc<- median(Test$YcenC)
   ZnP2<-function(A,Ncomp){
     Obj<- new("CmplxIm", img=A)
     momentType(Obj)<-"gpzm"
     setOrder(Obj)<-Ncomp
     setParams(Obj)<-1
     Moments(Obj)<- NULL
     Invariant(Obj) =NULL
     Reconstruct(Obj)<- c(dim(A)[1],dim(A)[2])
     return(Obj)
   }
   DimZ<-dim(Images)[3]

   Noi1<-array(dim=c(1,DimZ))
   Noi2<-Noi1
   Rout<-65
   Rin<-25
   Rout2<-35

   TF<-matrix("FALSE",nr=dim(Images)[1],nc=dim(Images)[2])
   xy<-expand.grid( (Xc-Rout):(Xc+Rout),(Yc-Rout+1):(Yc+Rout+1))
   dis<-sqrt( (xy[,1]-Xc)**2 + (xy[,2]-Yc)**2)
   wher<-which( dis>Rin & dis<Rout)
   xy.n<-xy[wher,]
   TF[xy.n[,1],xy.n[,2]]<-"TRUE"

   Ml1<-array(dim=c(Ncomp,Ncomp,DimZ))
    for(di in 1:DimZ){
     A0<-Images[(Xc-Rout2):(Xc+Rout2),(Yc-Rout2):(Yc+Rout2),di]
     A<-A0
     NasRe<-which(is.na(A)==TRUE | is.nan(A)==TRUE | abs(A)==Inf)
     if(length(NasRe)!=0){
       A[NasRe]<-0
     }
     k1<-ZnP2(A,Ncomp-1)
     Ml1[,,di]<-k1@moments
     kk0<-Images[,,di]
     nans<-which(is.na(kk0)==TRUE)
     if(length(nans)!=0){kk0[nans]<-0}
     Noi1[1,di]<-sd(kk0[TF=="TRUE"] )
    }
   Mr1<-matrix(0,nr=Ncomp,nc=DimZ)

   for(i in 1:Ncomp){
     for(k in 1:DimZ){
      if(i==1){
        Mr1[i,k]<-sqrt(sum(Mod(Ml1[,i,k])**2))/Noi1[1,k] ### @15/01/2018 I've changed this def. of modulus
      }else{
        Mr1[i,k]<-sqrt(sum(Mod(Ml1[,i,k])**2)) ### @15/01/2018 I've changed this def. of modulus
      }
   }}
   MST1<-matrix(nr=0,nc=Ncomp)
   Names.list<-matrix(nr=0,nc=1)
    for(j in 1:DimZ){
      MST1<-rbind(MST1,Mr1[,j])
    #  MST2<-rbind(MST2,Mr2[,j])
      Names.list<-rbind(Names.list,paste0(Name,"__",j))
    }

   MST1<-MST1
   MSQ1<-sqrt( MST1[,2]**2 + MST1[,3]**2 )
   MS01<-MST1[,1]

   MTF<-matrix(Names.list,nc=1)
   MTF<-cbind(MTF,MS01) # SNR_left
   MTF<-cbind(MTF,MSQ1) # Homog_left

   return(MTF)
  }


  EVA<-function(path.data,work.area){

    RT<-read.table("Super_F_table.txt",header=TRUE)
    SC1<-which(RT$Flag=="Scs")
    SK1<-which(RT$Flag=="Sky")

    names<-paste0("Scs_", gsub(".fits","_DF.fits",as.character(RT$FILE[SC1])) )

    Ncomp<-5
    A<-RT[SC1,]
    VU<-0
    Header<-c(names(RT),"SubName","Intensity","Homog")
    for(i in 1:length(names)){
      SM<-matrix(nr=0,nc=length(A[1,])+3)
      B<-readFITS(paste0(work.area,"/Data/",names[i]))$imDat
      M<-ZNPcNC(Ncomp,gsub('.fits','',names[i]),B)
      Et<-matrix(nr=0,nc=length(A[1,]))  
      for(tt in 1:dim(B)[3]){
        Et<-rbind(Et,A[i,])
      }
      SM<-rbind(SM,cbind(Et,M)  )
      rm(B)
      cat(paste0("------ ",i,"/",length(names),"------\n"))
      if(i==1){
        header<-matrix(Header,nr=1,nc=1)
        write.table(header,paste0(work.area,"/Total_inf.txt"),col.names=FALSE,row.names=FALSE,quote=FALSE,append=FALSE)
        write.table(SM,paste0(work.area,"/Total_Final_Info_NoRed.txt"),col.names=FALSE,row.names=FALSE,quote=FALSE,append=FALSE)
      }else{
       write.table(SM,paste0(work.area,"/Total_Final_Info_NoRed.txt"),col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE)
      }
    }
  }
   
######
######
if(EvA_File=="TRUE"){
 cat(" ---- Running the EVA file (Homog. and SNR Torus)...\n")
 EVA(path.data2,work.area2)
}
######
######

write.table("ok3","ok1.txt",col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE)


#/usr/lib64/R/bin/Rscript

#Nam<-c( "FILE","DPR.CATG","DPR.TYPE","MJD.OBS","EXPTIME","NAXIS1","NAXIS2","NEXPO","RA",          
# "DEC","Time","FWHM_start","FWHM_end","FWHM_IA","Coher_time","Wind_speed","Airmass_start",
#"Airmass_end","Focu_value","DimZ","NEXP","Dist","ObjType","NewName","Flag","numImage","XcenC","YcenC","SubName","Intensity","Homog")

Nam<-c( "FILE","DPR.CATG","DPR.TYPE","MJD.OBS","EXPTIME","NAXIS1","NAXIS2","NEXP","RA",          
 "DEC","Time","FWHM_start","FWHM_end","FWHM_IA","Coher_time","Wind_speed","Airmass_start",
"Airmass_end","Focu_value","DimZ","NEXPO","Dist","ObjType","NewName","Flag","numImage","XcenC","YcenC","SubName","Intensity","Homog")


FM<-read.table(paste0(work.area,"/Total_Final_Info_NoRed.txt"),header=FALSE)
names(FM)<-Nam


FM$Intensity<- (FM$Intensity - min(FM$Intensity))/(max(FM$Intensity)-min(FM$Intensity))

FM$Homog <- (-FM$Homog  + max(FM$Homog) )/(max(FM$Homog) - min(FM$Homog)  )


write.table(FM,paste0(work.area,"/Total_Final_Info_Red.txt"),col.names=TRUE,row.names=FALSE,quote=FALSE)













