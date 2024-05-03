#find volume of a backwater, used to simplify calculations for segments that don't fully contain the BW (closely spaced jams)
findBWvolv2<-function(Lsi,Bi,colh0,colhji,HBFi,Beffi,Si,Qi) {

  colVchannel <- (colhji>colh0)*findViterate(Lsi,Bi,colh0,colhji,HBFi,Beffi,Si,Qi)
  
  matMin=matrix(,length(colh0),2)
  matMin[,1]=colhji-HBFi
  matMin[,2]=colhji-colh0
  minDeltaHFP=(colhji>HBFi)*rowMins(matMin,TRUE)
  matMin=matrix(,length(colh0),2)
  matMin[,1]=colh0-HBFi
  matMin[,2]=0
  minDeltaH0FP=(colhji>HBFi)*rowMaxs(matMin,TRUE)
  
  matMin=matrix(,length(colh0),2)
  matMin[,1]=HBFi
  matMin[,2]=colh0
  whatisHFP=(colhji>HBFi)*rowMaxs(matMin,TRUE)
  
  colVFP <- (colhji>HBFi)*((Beffi-1)*Bi/(2*Si)*(minDeltaHFP)^2) #iterate
  colVFP[is.nan(colVFP)] = 0
  
  BWvol <- colVchannel+colVFP
  return(BWvol)
  
}