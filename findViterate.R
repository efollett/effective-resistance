findViterate<-function(Lsi,Bi,colh0,colhji,HBFi,Beffi,Si,Qi) {
  #calculates asymptotic backwater curve
  allpts=seq(1,length(colhji),1)
  colVI=matrix(,nrow = 1,ncol = length(colhji))
  
  for (p in allpts){
   HbwX=colhji[p]
   delX=0.1
   hC=((Qi[p]/Bi)^2/9.8)^(1/3)
   LbwX=0
   VbwX=0
   #standard backwater case
   while (HbwX>colh0[p]*1.01) {
    HbwX1=HbwX-delX*Si*(1-(colh0[p]/HbwX)^3)/(1-(hC/HbwX)^3)
    VbwX=VbwX+delX*0.5*(HbwX1+HbwX)*Bi
    LbwX=LbwX+delX
    HbwX=HbwX1
   }
   #overlapping backwater case
   if (LbwX>Lsi) {
     LbwX=0
     VbwX=0
     HbwX=colhji[p]
     while  (LbwX < (Lsi-delX)) {
       HbwX1=HbwX-delX*Si*(1-(colh0[p]/HbwX)^3)/(1-(hC/HbwX)^3)
       VbwX=VbwX+delX*0.5*(HbwX1+HbwX)*Bi
       LbwX=LbwX+delX
       HbwX=HbwX1

     }
   }
  colVI[1,p]=(VbwX-colh0[p]*Bi*LbwX) #report excess volume above uniform flow only
  }
  return(colVI)
}