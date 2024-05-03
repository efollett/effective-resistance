#find volume of uniform flow in segment, used to simplify for calculation of closely spaced jams
findunifvol <- function(Lsi,Bi,colh0,colhj,HBFi,Beffi) {
  
  colVchannel <- Lsi*Bi*colh0

  matMin=matrix(,length(colh0),2)
  matMin[,1]=colhj-HBFi
  matMin[,2]=colhj-colh0
  minDeltaHFP=(colhj>HBFi)*rowMins(matMin,TRUE)
  
  matMin=matrix(,length(colh0),2)
  matMin[,1]=colh0-HBFi
  matMin[,2]=0
  minDeltaH0FP=(colhj>HBFi)*rowMaxs(matMin,TRUE)

  colVFP <- (colhj>HBFi)*((Beffi-1)*Bi*Lsi*minDeltaH0FP)
  colVFP[is.nan(colVFP)] = 0
  
  unifvol <- colVchannel+colVFP
  return(unifvol)
}