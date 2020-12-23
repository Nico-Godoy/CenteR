
memory.limit(32768)

work.area<-work.area2
cl <- makeCluster(Cluster_number)
registerDoParallel(cl)


Object<-"Test"
work.area<-work.area2
pdc<-path.data2
pdp<-paste0(work.area2,"Data")

write.table(date(),paste0(work.area,"Date_prep_input.txt"),col.names=FALSE,row.names=FALSE,quote=FALSE,append=FALSE)


if(HDF5_file_tot=="TRUE"){

Transf<-function(work.area){ # Estimado la distancia entre AGPM y Circle
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
   return(list(p=c(XT,YT),err=c(DXT,DYT)))
}

Trs<-Transf(work.area)
XT<- Trs$p[1]  # -10.58339  10.38168
YT<- Trs$p[2]

List<-read.table(paste0(work.area,"Final_table_complete_information_v6.txt"),
                 header=TRUE)


Dist<-sqrt((List$Sy - List$CorY)**2 + (List$Sx - List$CorX)**2  )
sel<-which(Dist<=median(Dist)+3*mad(Dist) &  
           List$RES<=median(List$RES)+3*mad(List$RES) & List$Bkg>=0  &
           List$index!=1 & List$E_SDx/List$SDx<re & List$E_SDy/List$SDy<re & 
           List$E_Rho/List$Rho<re & List$E_CoS/List$CoS<re &
           List$BkgSD<= median(List$BkgSD) +3*mad(List$BkgSD)  #) #&
          & List$Int>median(List$Int) + SiG*mad(List$Int))

sel<-1:length(List[,1])


NoData<-which(List$Xc>=99999 | List$Yc>=99999 | List$Sx==9999)
if(length(NoData)!=0){
   sel<-sel[-NoData]
   cat(paste0("The bad data correspond to: ",length(NoData),"\n"))
}

#################
Names_1<-unique(List$FILE[sel])

val_1<-0
for(na1 in 1:length(Names_1)){
  sel_1<-which(List$FILE[sel]==Names_1[na1])
  aux_1<-which(List$index[sel[sel_1]]==List$DimZ[sel[sel_1]][1])
  if(length(aux_1)!=0){
   val_1[na1]<-sel_1[aux_1]
  }else{
   val_1[na1]<-0
  }
}

if(length(which(val_1==0))!=0){
  val_1<-val_1[-which(val_1==0)]
}

if(length(val_1)!=0){
   sel<-sel[-val_1]
}

#################


# ~14% of data
cat(paste0("---",100*length(sel)/dim(List)[1]," % of data \n"  ))
cat("All is ok @ this point 1\n")

################
# Added to use the same values as in the paper!!
# ValAux<-read.table("/home/ngodoy/data/NACO_data/Processing/BetaPic_p/Vals_match.ascii")
# sel<-sel[ValAux[,1]]

################

write.table(c(100*length(sel)/dim(List)[1],length(sel),length(List[,1])),paste0(work.area,"Perc_co.txt"),col.names=FALSE,row.names=FALSE)


subDir<-strsplit(work.area2,"/")[[1]]

work.area3<-paste0("/debris_data/ngodoy/PynPoint_products/",subDir[length(subDir)])

system(paste0("mkdir ",work.area3))


name.obj<-paste0(work.area3,"/",Object_hdf5_tot)
File<-h5file(name=name.obj,mode="a")
#File["/header_im_arr/NEW_PARA"]<-matrix(nr=1,nc=length(sel))
File["/header_im_arr/Used_Files"]<-matrix(nr=1,nc=length(sel))
#File["/im_arr"]<-array(dim=c(Xl,Yl,length(sel)))
createGroup(File,"/Object_V5.hdf5")
createGroup(File,"/hesder_im_arr")
createGroup(File,"/hesder_im_arr/NEW_PARA")

cat("All is ok @ this point 2\n")

Names<-unique(List$FILE[sel])


if(Re_Size=="TRUE"){

###########
X_mean<- mean(List$Sx[sel],na.rm=T)
Y_mean<- mean(List$Sy[sel],na.rm=T)

     if( abs( X_mean-floor(X_mean) - 0.5 )<=0.25 ){
          if(Xl%%2 !=0){ 
            Xl<-Xl+1
          } 
     }else{
          if(Xl%%2 !=1){ 
            Xl<-Xl+1
          }
     }   
     if( abs( Y_mean-floor(Y_mean) - 0.5 )<=0.25 ){
          if(Yl%%2 !=0){ 
            Yl<-Yl+1
          }
     }else{
          if(Yl%%2 !=1){ 
            Yl<-Yl+1
          }
     }

# Cuadratura
     if(Xl%%2==0){
        Xl<-Xl+1
     }

     if(Yl!=Xl){
       if(  abs(X_mean - round(X_mean)) > abs(Y_mean - round(Y_mean))  ){
          Xl<-Yl
       }else{
          Yl<-Xl
       }
     }

###########


}
###########

Xl<-Xl+4
Yl<-Yl+4


File["/im_arr"]<-aperm(foreach(i=1:length(Names), .combine=abind,.packages=c('foreach','abind','FITSio','fields','akima'),.export=c("pdc","pdp","Xl","Yl"))%dopar%{
     s2<-which(as.character(Names[i])==List$FILE[sel])
     Num<-List$index[sel[s2]]
     Name0<-strsplit(as.character(List$FILE[[sel[s2[1]]]]),".fits")[[1]]
     A<-readFITS(paste0(pdp,"/Scs_",Name0,"_DF.fits"))$imDat
     A0<-array(dim=c(Xl,Yl,length(Num)))
     for(j in 1:length(Num)){
#       Sx<-List$Sx[sel[s2[j]]]
#       Sy<-List$Sy[sel[s2[j]]]
       Sx<- List$CorX[sel[s2[j]]]  #List$Sx[sel[s2[j]]]
       Sy<- List$CorY[sel[s2[j]]]  #List$Sy[sel[s2[j]]]
###### New XR and YR
#       if( abs( X_mean-floor(X_mean) - 0.5 )<=0.25 ){
#          XR<- seq( floor(Sx) - floor(Xl/2) +1,floor(Sx) + floor(Xl/2) ,by=1 )
#       }else{
#          XR<- ( ( round(Sx) - floor(Xl/2)):(round(Sx) + floor(Xl/2)) ) 
#       }   
#       if( abs( Y_mean-floor(Y_mean) - 0.5 )<=0.25 ){
#          YR<- seq( floor(Sy) - floor(Yl/2) +1,floor(Sy) + floor(Yl/2) ,by=1 ) 
#       }else{
#          YR<- ( round(Sy) - floor(Yl/2)):(round(Sy) + floor(Yl/2))  
#       }
######

       XR<- seq( round(Sx) - floor(Xl/2) -5,round(Sx) + floor(Xl/2) +5,by=1 )
       YR<- seq( round(Sy) - floor(Yl/2) -5,round(Sy) + floor(Yl/2) +5,by=1 )

       XRn<- seq( Sx - floor(Xl/2) , Sx + floor(Xl/2) ,by=1 )
       YRn<- seq( Sy - floor(Yl/2) , Sy + floor(Yl/2) ,by=1 )

#       AuxIm<-image.smooth(A[XR,YR,Num[j]], xwidth=1.25,ywidth=1.25)$z
       AuxIm<-image.smooth(A[XR,YR,Num[j]],xwidth=1.25,ywidth=1.25)$z

       AuxImN<-bicubic(x=XR,y=YR,z=AuxIm,x0=expand.grid(XRn,YRn)[,1] , y0=expand.grid(XRn,YRn)[,2] )

       Aux01<- matrix(AuxImN$z,nc=length(XRn))

       A0[,,j]<-t(Aux01)
     }
     NasRem<-which(abs(A0)==Inf | is.na(A0)==TRUE | is.nan(A0)==TRUE)
     if(length(NasRem)!=0){
       A0[NasRem]<-median(A0[-NasRem])
     }
     return(A0)
})

File["/header_im_arr/Shift"]<-foreach(i=1:length(Names), .combine=rbind,.packages=c('foreach','abind','FITSio','akima','fields'),.export=c("pdc","pdp","Xl","Yl","work.area"))%dopar%{
     s2<-which(as.character(Names[i])==List$FILE[sel])

       SCx<- List$Sx[sel[s2]] - List$CorX[sel[s2]]
       SCy<- List$Sy[sel[s2]] - List$CorY[sel[s2]]
 
       SCxy<-matrix(c(SCx,SCy),nc=2)

     return(SCxy)
}



cat("All is ok @ this point 3\n")

#File["/im_arr"]<-aperm(Aux2)
#rm(Aux2)

cat("All is ok @ this point 4\n")


if(PAng=="Old"){
Aux<-foreach(i=1:length(Names), .combine=cbind,.packages=c('foreach','FITSio','astro'),.export=c("pdc","pdp","Xl","Yl"))%dopar%{
   s2<-which(as.character(Names[i])==List$FILE[sel])
   Num<-List$index[sel[s2]]
     Name0<-strsplit(as.character(Names[i]),".fits")[[1]]
     hdr<-read.fitshdr(paste0(pdc,"/",Name0,".fits"))
     DIT<-as.numeric(hdr[hdr[,1]=="ESO DET DIT",2])
     NDIT<-as.numeric(hdr[hdr[,1]=="ESO DET NDIT",2])
     DIMZ<-dim(readFITS(paste0(pdc,"/",Name0,".fits"))$imDat)[3] -1
  if(is.na(DIMZ)!=TRUE & length(NDIT)!=0){
   if(DIMZ!=1){

     Ini<-as.numeric(hdr[hdr[,1]=="ESO ADA POSANG",2])
     End<-as.numeric(hdr[hdr[,1]=="ESO ADA POSANG END",2])
     if(NDIT>DIMZ){ ZL<-DIMZ}else{ ZL<-NDIT}
     PA<-seq(Ini,End,len=ZL)

     Name<-Names[i]
     aux0<-strsplit(as.character(Name),".fits")
     aux1<-strsplit(aux0[1][[1]][1],"T")
     aux2<-as.numeric(strsplit(aux1[1][[1]][2],":")[1][[1]])
     time<- aux2[1] + aux2[2]/60 + aux2[3]/3600
     Time<- seq( time, time + ZL*DIT/3600,by=DIT/3600)[-c(ZL+1)]

     Para<-matrix(-PA[Num],nr=1)
     return(Para)
  }}
}
}else{
library(rPython)
Aux<-foreach(i=1:length(Names),.combine=cbind,.packages=c('foreach','FITSio','astro','rPython','astroFns'),.export=c("pdc","pdp","Xl","Yl"))%dopar%{
     s2<-which(as.character(Names[i])==List$FILE[sel])
     Num<-List$index[sel[s2]]
     Name0<-strsplit(as.character(List$FILE[[sel[s2[1]]]]),".fits")[[1]]
   DIMZ<-dim(readFITS(paste0(pdp,"/Scs_",Name0,"_DF.fits"))$imDat)[3] -1
     hdr<-read.fitshdr(paste0(pdc,"/",Name0,".fits"))
     new_angles<-matrix(nr=0,nc=1)
     pupil_pos_arr<-matrix(nr=0,nc=1)
     # pupil offset in degrees
     m_pupil_offset <- 0.            # No offset here
     # no overheads in cube mode, since cube is read out after all individual exposures
     # see NACO manual page 62 (v102)
     m_O_START <- 0.
     m_DIT_DELAY <- 0.
     m_ROT<- 0.
     # rotator offset in degrees
     m_rot_offset = 89.44 # According to NACO manual page 65 (v102)

     ndit<- min( c(as.numeric(hdr[hdr[,1]=="ESO DET NDIT",2]), DIMZ )) # "NFRAMES" HIERARCH ESO DET NDIT
     exptime<- as.numeric(hdr[hdr[,1]=="ESO DET DIT",2])/3600 #
     tel_lat<- as.numeric(hdr[hdr[,1]=="ESO TEL GEOLAT",2])
     tel_lon<- as.numeric(hdr[hdr[,1]=="ESO TEL GEOLON",2])

#     ra<- as.numeric(hdr[hdr[,1]=="RA",2]) # using all, compute the mean
#     dec<- as.numeric(hdr[hdr[,1]=="DEC",2]) # using all, compute the mean
     ra<-mean(List$RA,na.rm=T)
     dec<-mean(List$DEC,na.rm=T)
     obs_dates<- hdr[hdr[,1]=="DATE",2]
     pupil_pos<- as.numeric(hdr[hdr[,1]=="ESO ADA PUPILPOS",2])

     python.exec("import ephem")
     python.exec("obs = ephem.Observer()")
     python.exec(paste0("obs.lat = ephem.degrees(str('",tel_lat,"'))"))
     python.exec(paste0("obs.long = ephem.degrees(str('",tel_lon,"'))"))
     python.exec(paste0("obs.date = '",as.character(gsub("T"," ",obs_dates[[1]])),"'"))
     python.exec("print obs.date")
     # Get sideral time in hours
     python.exec("sid_time = str(obs.sidereal_time())")
     # Get hours minutes and seconds
     python.exec("h, m, s = sid_time.split(':')")
     python.exec("sid_time_2 = (float(h) + (float(m) / 60.) + (float(s) / 3600.))")
     sid_time_py<-python.get("sid_time_2")
     sid_time_arr <- seq(sid_time_py+m_O_START,
            (sid_time_py + m_O_START) + (exptime + m_DIT_DELAY + m_ROT)*(ndit),
             len=ndit)
     # Convert to degrees
     sid_time_arr_deg <- sid_time_arr * 15.
     # Calculate hour angle in degrees
     hour_angle <- sid_time_arr_deg - ra
     # Conversion to radians:
     hour_angle_rad <- dms2rad(hour_angle)
     dec_rad <- dms2rad(dec)
     lat_rad <- dms2rad(tel_lat)
     p_angle <- atan2( sin(hour_angle_rad),
                 (cos(dec_rad)*tan(lat_rad)-sin(dec_rad)*cos(hour_angle_rad)))
     new_angles <- rbind(new_angles, matrix(p_angle*180/pi,nc=1))
     pupil_pos_arr <- rbind(pupil_pos_arr, matrix(pupil_pos,nr=ndit,nc=1))
     # Correct for pupil offset (NACO)
     # See NACO manual page 65 (v102)
     new_angles_corr <- new_angles - (90. + (m_rot_offset-pupil_pos_arr))
     PARAng<-new_angles_corr[,1]
     Para<-matrix(-PARAng[Num],nr=1)
     return(Para)
}

}
cat("All is ok @ this point 5\n")

File["/header_im_arr/RESP_PARANG"]<-as.vector(Aux[1,])
File["/header_im_arr/PARANG"]<-as.vector(Aux[1,])*0
rm(Aux)

cat("--------------------\n")
cat(paste0("dim Im: ",dim(File["/im_arr"][]),"  --  len PA: ",length(File["/header_im_arr/PARANG"][]),"\n"))
cat("--------------------\n")


cat("All is ok @ this point 6\n")
createAttribute(File["/im_arr"], "PIXSCALE", 0.02719)
h5close(File)

stopCluster(cl)

system(paste0("cp ",name.obj," ",work.area))


write.table(date(),paste0(work.area,"Date_prep_input.txt"),col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE)

}else{
   Xl<-Xl+4
   Yl<-Yl+4
}


############################################################################################################################################
############################################################################################################################################


if(PCA_run_tot=="TRUE"){

system(paste0("mkdir ",work.area2,HDF5_folder_tot))
system(paste0("mkdir ",work.area2,HDF5_folder_tot,"/Products"))
system(paste0("cp ",work.area,Object_hdf5_tot," ",work.area2,HDF5_folder_tot,"/Object_test.hdf5"))

subDir<-strsplit(work.area2,"/")[[1]]

work.area3<-paste0("/debris_data/ngodoy/PynPoint_products/",subDir[length(subDir)],"/")

system(paste0("mkdir ",work.area3))
system(paste0("mkdir ",work.area3,HDF5_folder_tot))
system(paste0("mkdir ",work.area3,HDF5_folder_tot,"/Products"))
system(paste0("cp ",work.area,Object_hdf5_tot," ",work.area3,HDF5_folder_tot,"/Object_test.hdf5"))


cat("Running multiple PCA ...\n")

memory.limit(32768)


Path1Zero<-paste0(work.area2,HDF5_folder_tot,"/Products")
Path1<-paste0(work.area3,HDF5_folder_tot,"/Products")
Path2<-paste0(work.area3,HDF5_folder_tot)

#COMP<-unique(round(10**seq(log10(1),log10(23400*0.5),len=30)))
Radius<-Mask_radius

#COMP<-c(65,85,95,110)

  system(paste0("mkdir ",Path1))


for(jj in 1:length(COMP)){
  Path3<-paste0(work.area3,HDF5_folder_tot,"/Temporal_Folder_",jj)
  Path4<-paste0(Path3,"/Products")

  system(paste0("mkdir ",Path3))
  system(paste0("mkdir ",Path4))
  system(paste0("rm ",Path4,"/PynPoint_database.hdf5"))
  system(paste0("cp ",path.pip,"PynPoint_processing__A.py ",Path3,"/PP.py"))
  Object<-paste0(Path4,"/Image_",COMP[jj],"_PC.png")
  system(paste0("sed -i \'s|PATH1|\"",Path4,"\"|gi\' ",Path3,"/PP.py"))
  system(paste0("sed -i \'s|PATH2|\"",Path2,"\"|gi\' ",Path3,"/PP.py"))
  system(paste0("sed -i \'s|COMP|",COMP[jj],"|gi\' ",Path3,"/PP.py"))
  system(paste0("sed -i \'s|RADIUS|",Radius,"|gi\' ",Path3,"/PP.py"))
  system(paste0("sed -i \'s|OBJECT|\"",Object,"\"|gi\' ",Path3,"/PP.py"))
  system(paste0("chmod +x ",Path3,"/PP.py"))

  Path3C<-paste0(work.area3,HDF5_folder_tot,"/SUB_Temp_Folder_",jj)
  Path4C<-paste0(Path3C,"/Products")
  Path2C<-paste0(Path3C)

  system(paste0("mkdir ",Path3C))
  system(paste0("mkdir ",Path4C))
  system(paste0("rm ",Path4C,"/PynPoint_database.hdf5"))
  system(paste0("cp ",path.pip,"PynPoint_processing__B.py ",Path3C,"/PP.py"))
  system(paste0("sed -i \'s|PATH1|\"",Path4C,"\"|gi\' ",Path3C,"/PP.py"))
  system(paste0("sed -i \'s|PATH2|\"",Path2C,"\"|gi\' ",Path3C,"/PP.py"))

  system(paste0("chmod +x ",Path3C,"/PP.py"))

}

cl <- makeCluster(Cluster_number_p2)
registerDoParallel(cl)

KK<-foreach(jjj=1:length(COMP), .combine=rbind,.packages=c('foreach','abind','FITSio','fields','akima','h5','astro'),.export=c("Radius","Path1","Path2","COMP","Object"),.inorder=TRUE)%dopar%{

 ZeroFn<-function(COMPT){
   Len0<-length(strsplit(as.character(COMPT), split="")[[1]])
   Val<-6-Len0
   Zeros<-paste0(as.character(seq(0,0,len=Val)),collapse="") 
   return(Zeros)
 }

 if(file.exists(paste0(Path1,"/Image_",ZeroFn(COMP[jjj]),COMP[jjj],"_PC_rz.fits"))==TRUE  ){
   FileLEN<-as.numeric(strsplit(  system(paste0("ls -lF ",Path1,"/Image_",ZeroFn(COMP[jjj]),COMP[jjj],"_PC_rz.fits"),intern=TRUE), split=" ")[[1]][5])
 }else{
   FileLEN<-0
 }

# if(file.exists(paste0(Path1,"/Image_",ZeroFn(COMP[jjj]),COMP[jjj],"_PC_rz.fits"))!=TRUE | FileLEN==0 ){

    Path3<-paste0(work.area3,HDF5_folder_tot,"/Temporal_Folder_",jjj)
    Path4<-paste0(Path3,"/Products")

    system(paste0("python ",Path3,"/PP.py"))

}

stopCluster(cl)


cl <- makeCluster(Cluster_number_p2)
registerDoParallel(cl)


KK<-foreach(jjj=1:length(COMP), .combine=rbind,.packages=c('foreach','abind','FITSio','fields','akima','h5','astro'),.export=c("Radius","Path1","Path2","COMP","Object"),.inorder=TRUE)%dopar%{

 ZeroFn<-function(COMPT){
   Len0<-length(strsplit(as.character(COMPT), split="")[[1]])
   Val<-6-Len0
   Zeros<-paste0(as.character(seq(0,0,len=Val)),collapse="") 
   return(Zeros)
 }

    Path3<-paste0(work.area3,HDF5_folder_tot,"/Temporal_Folder_",jjj)
    Path4<-paste0(Path3,"/Products")

    A<-h5file(paste0(Path4,"/PynPoint_database.hdf5"),mode="r")

    PSF_sub<-aperm(A[paste0("/PSF_sub",COMP[jjj])][]) # array of  Lx x Ly x numIm
    PA<-A["/header_im_arr/RESP_PARANG"][] # vector of L_numIm
#    Shift<- h5file(paste0(work.area3,HDF5_folder,"/Object_test.hdf5"),mode="r")["/im_arr_shift"][] # matrix of numIm x [x,y]
    Shift<- A["/header_im_arr/Shift"][] # matrix of numIm x [x,y]

    Xl0<-dim(PSF_sub)[1]
    Yl0<-dim(PSF_sub)[2]

    PSF_sub_shift<-array(dim=c( dim(PSF_sub)[1]-4 ,dim(PSF_sub)[2]-4,dim(PSF_sub)[3]   ))

library(akima)
library(fields)

    for(tr in 1:dim(PSF_sub)[3]){
       XR<- seq( - floor(Xl0/2) ,  floor(Xl0/2) ,by=1 )
       YR<- seq( - floor(Yl0/2) ,  floor(Yl0/2) ,by=1 )

       XRn<- seq( Shift[tr,2]*0 - floor((Xl0-4)/2) , Shift[tr,2]*0 + floor((Xl0-4)/2) ,by=1 ) # Ive changed the 2 by 1 in shift index , NO # becuase Ive saved the transp of the image! 
       YRn<- seq( Shift[tr,1]*0 - floor((Yl0-4)/2) , Shift[tr,1]*0 + floor((Yl0-4)/2) ,by=1 )

       AuxImN<-bicubic(x=XR,y=YR,z=PSF_sub[,,tr],x0=expand.grid(XRn,YRn)[,1] , y0=expand.grid(XRn,YRn)[,2] )
       Aux01<- matrix(AuxImN$z,nc=length(XRn))
     
       PSF_sub_shift[,,tr]<- Aux01

    }

  Path3C<-paste0(work.area3,HDF5_folder_tot,"/SUB_Temp_Folder_",jjj)
  Path4C<-paste0(Path3C,"/Products")
  Path2C<-paste0(Path3C)

   system(paste0("mkdir ",Path3C))

    name.objC<-paste0(Path3C,"/Object_PSF_sub_shifted.hdf5")
  system(paste0("rm ",name.objC))
    FileC<-h5file(name=name.objC,mode="a")
    FileC["/for_derotate"]<-aperm(PSF_sub_shift)
    FileC["/header_for_derotate/PARANG"]<-PA
    createAttribute(FileC["/for_derotate"], "PIXSCALE", 0.02719)
    h5close(FileC)

}


stopCluster(cl)


cl <- makeCluster(Cluster_number_p2)
registerDoParallel(cl)


KK<-foreach(jjj=1:length(COMP), .combine=rbind,.packages=c('foreach','abind','FITSio','fields','akima','h5','astro'),.export=c("Radius","Path1","Path2","COMP","Object"),.inorder=TRUE)%dopar%{

 ZeroFn<-function(COMPT){
   Len0<-length(strsplit(as.character(COMPT), split="")[[1]])
   Val<-6-Len0
   Zeros<-paste0(as.character(seq(0,0,len=Val)),collapse="") 
   return(Zeros)
 }

  Path3C<-paste0(work.area3,HDF5_folder_tot,"/SUB_Temp_Folder_",jjj)
  Path4C<-paste0(Path3C,"/Products")
  Path2C<-paste0(Path3C)

  system(paste0("rm ",Path4C,"/PynPoint_database.hdf5"))
  system(paste0("rm ",Path3C,"Products/*"))
   system(paste0("python ",Path3C,"/PP.py"))

   AFS<-h5file(paste0(Path4C,"/PynPoint_database.hdf5"),mode="r")
   Der_im <- AFS["/derotated_im"][]

   writeFITSim(Der_im,paste0(Path1Zero,"/Image_",ZeroFn(COMP[jjj]),COMP[jjj],"_PC_rz_nr_test.fits"))

}

stopCluster(cl)

################## smoothing


if(1==2){

cat("Smooting the FITS-PCA ...\n")

library(fields)
library(FITSio)

path<-Path1Zero

system(paste0("ls ",path,"/*.fits > ",path,"/Files.list"))
system(paste0("sed -i \'s|",path,"||g\' ",path,"/Files.list"))

List<-read.table(paste0(path,"/Files.list"))

system(paste0("mkdir ",path,"/Smooth_files"))
for(i in 1:length(List[,1])){
  Nam<-as.character(List[i,1])
  IM<-readFITS(paste0(path,"/",Nam))$imDat
  IMs<- image.smooth(IM, xwidth=1.25,ywidth=1.25)
  Name<-gsub(".fits","",Nam)
  writeFITSim(IMs$z,paste0(path,"/Smooth_files/",Name,"_smoothed.fits"))
#  png(paste0(path,"/Smooth_files/",Name,"_smoothed.png"))
#   par(oma=c(0.8, 0, 0, 0.8))
#   image.plot(IMs)
#  dev.off()
}
}

}


