
UNIFIER<-function(work.area){
   List<-read.table(paste0(work.area,"Final_table_complete_information_v3.txt"),header=TRUE)
   L2<-read.table(paste0(work.area,"Total_Final_Info_NoRed.txt"))
   D<-length(L2[1,])
   List2<-matrix(nr=0,nc=2)
   for(i in 1:length(List[,1])){
     aux<-which(List$SubName[i]==L2[,D-2])
     List2<-rbind(List2,L2[aux,D-c(1,0)])
   }
   nams<-colnames(List)
   Nams<-c(nams,"ZP_Inten","ZP_Homog")

   List3<-cbind(List,List2)
    colnames(List3)<-Nams
   write.table(List3,paste0(work.area,"Final_table_complete_information_v3_ZP.txt"),col.names=TRUE,row.names=FALSE,quote=FALSE)
}


