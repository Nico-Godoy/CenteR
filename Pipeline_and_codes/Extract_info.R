################################################################################

memory.limit(32768)


################################################################################


memory.limit(32768)

work.area<-work.area2
cl <- makeCluster(Cluster_number)
registerDoParallel(cl)


Object<-"Test"
work.area<-work.area2
pdc<-path.data2
pdp<-paste0(work.area2,"Data")




##
XT<- -10.45105
YT<- +11.71069

List<-read.table(paste0(work.area,"Final_table_complete_information_v3.txt"),
                 header=TRUE)

if(Frame_selec=="Yes"){
Dist<-sqrt((List$Sy - List$CorY)**2 + (List$Sx - List$CorX)**2  )
sel<-which(Dist<=median(Dist)+3*mad(Dist) &  
           List$RES<=median(List$RES)+3*mad(List$RES) & List$Bkg>=0  &
           List$index!=1 & List$E_SDx/List$SDx<re & List$E_SDy/List$SDy<re & 
           List$E_Rho/List$Rho<re & List$E_CoS/List$CoS<re &
           List$BkgSD<= median(List$BkgSD) +3*mad(List$BkgSD)  #) #&
          & List$Int>median(List$Int) + SiG*mad(List$Int) & List$DimZ>=100)
}else{
sel<-which(List$DimZ>100)
}

 X<-log10(List$RES[sel])
 Y<-log10(List$Nor1[sel]) 

   X0<-seq(min(X),max(X),len=100)
   Y1<-0
   X1<-0
   DY1<-0
   DX1<-0
   for(i in 1:(length(X0)-1)){
     aux<-which( X>=X0[i] & X<=X0[i+1])
     Y1[i]<- median(Y[aux])
     X1[i]<-median(X[aux])
     DY1[i]<-mad(Y[aux])
     DX1[i]<-mad(X[aux])
   }

NAs<-which(is.na(Y1)==TRUE | is.na(X1)==TRUE | abs(Y1)==Inf)
if(length(NAs)!=0){
  X1<-X1[-NAs]
  Y1<-Y1[-NAs]
  DX1<-DX1[-NAs]
  DY1<-DY1[-NAs]
}

Ys<-spline(X1,Y1,xout=X)$y
Yn<-Y-Ys
Xn<-X

Non<-which(Yn>=median(Yn)+7*mad(Yn))

if(length(Non)!=0){
  sel<-sel[-Non]
}
##


Dist<-seq(Dis0,(Xl-1)/2 -Dis0,by=18)

#write.table(c("LOG"),paste0(work.area,"/log_ing.ascii"),col.names=FALSE,row.names=FALSE,quote=FALSE)

Names<-unique(List$FILE[sel])

File<-foreach(iii=1:length(Names), .combine=rbind,.packages=c('foreach','abind','FITSio','fields','akima','h5','astro'),.export=c("Dist","PA0"),.inorder=TRUE)%dopar%{
###
     s2<-which(as.character(Names[iii])==List$FILE[sel])
     Num<-List$index[sel[s2]]
     Name0<-strsplit(as.character(List$FILE[[sel[s2[1]]]]),".fits")[[1]]
     #A<-readFITS(paste0(pdp,"/Scs_",Name0,"_DF.fits"))$imDat
     #A00<-array(dim=c(Xl,Yl,length(Num)))

     #NasRem<-which(abs(A)==Inf | is.na(A)==TRUE | is.nan(A)==TRUE)
     #if(length(NasRem)!=0){
     #  A[NasRem]<-median(A[-NasRem])
     #}

     hdr<-read.fitshdr(paste0(pdc,"/",Name0,".fits"))
     Ini<-as.numeric(hdr[hdr[,1]=="ESO ADA POSANG",2])
     End<-as.numeric(hdr[hdr[,1]=="ESO ADA POSANG END",2])
     PA<-seq(Ini,End,len=100)
     Para<-matrix(-PA,nr=1) # starting with frame #2

     Time<-min(List$Time[sel[s2]])
     Dtime<-min(List$EXPTIME[sel[s2]])
   
     SM<-matrix(nr=0,nc=6+1)


     for(j in 1:length(Num)){
       Sx<-List$Sx[sel[s2[j]]]
       Sy<-List$Sy[sel[s2[j]]]
       XR<-(round(Sx)-floor(Xl/2)):(round(Sx)+floor(Xl/2))
       YR<-(round(Sy)-floor(Yl/2)):(round(Sy)+floor(Yl/2))
       #Aux01<-A[XR,YR,Num[j]]

       #IM<-t(Aux01)
       #A0<-IM
       PAg<- -Para[1,]
       #X<-1:dim(IM)[1]
       #Y<-1:dim(IM)[2]
       for(ij in 1:length(Dist)){ #1:length(Dist)){
        # Xll<-1:dim(Photom)[1] - (dim(Photom)[1]-1)/2 -1
        # Yll<-1:dim(Photom)[2] - (dim(Photom)[2]-1)/2 -1

         Xi<- Dist[ij]*cos((PAg[Num[j]]+PA0)*pi/180) 
         Yi<- Dist[ij]*sin((PAg[Num[j]]+PA0)*pi/180)

        # Xir<- Xi - floor(Xi)
        # Yir<- Yi - floor(Yi) 

        # Xli<-Xll + Xir
        # Yli<-Yll + Yir

        # Xl0<-2:(dim(Photom)[1]-1) - (dim(Photom)[1]-1)/2 -1
        # Yl0<-2:(dim(Photom)[2]-1) - (dim(Photom)[2]-1)/2 -1
        # Photom<-(Photom/max(Photom))*Intensidad
       
        # Zi<-bilinear.grid(x=Xli,y=Yli,z=Photom,xlim=range(Xl0),ylim=range(Yl0),dx=1,dy=1)$z
        # mxr<-mean(X)+(floor(Xi)+Xl0[1]):(floor(Xi)+Xl0[length(Xl0)])
        # myr<-mean(Y)+(floor(Yi)+Yl0[1]):(floor(Yi)+Yl0[length(Yl0)])
        # A0[mxr,myr]<-A0[mxr,myr]+Zi
 
         vector0<-c(Name0,List$index[sel[s2[j]]],Time+Dtime*List$index[sel[s2[j]]]/100,PAg[Num[j]],Xi,Yi,Dist[ij])
         SM<-rbind(SM,vector0)

       }
       #A00[,,j]<-A0

     }


#     write.table(paste0("--",i,"/",length(Names),"--"),paste0(work.area,"/log_ing.ascii"),col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE)
     return(SM)
}

#####


cat("All is ok @ this point 3\n")

#File["/im_arr"]<-aperm(Aux2)
#rm(Aux2)


colnames(File)<-c("Name","index","Time","PA","X","Y","Dist")

write.table(File,paste0(work.area,"/Summary_info.txt"),col.names=TRUE,row.names=FALSE,quote=FALSE)

cat("All is ok @ this point 4\n")


stopCluster(cl)





############################################################################################################################################
############################################################################################################################################



