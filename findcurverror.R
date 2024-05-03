
#find sum of residuals^2 to minimise with fminsearch
findcurverror <- function(coeffGuess,alloutJQ,riv){
  
  #find outflow curve WITH NO JAMS (isjam==0)
  alloutN <- findoutflowv7(riv=riv,coeff=coeffGuess,jams=0)
  residuals=alloutJQ-alloutN$Q #choose
  error=sum(residuals^2)
  disp(coeffGuess,error)
return(error)
}