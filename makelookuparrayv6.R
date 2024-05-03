makelookuparrayv6 <- function(CA,S,coeff,B,Ls,N,J,a,Hj,HBF,Beff,QBF,nFP) {
  library(Rfast)
  #---
  Hpts=seq(0,10,0.01) #Range of values for HJ/HBF--limit is set unreasonably high so that values always available
  rivarray <- array(0,dim=c(length(Hpts),5,N)) 
  Npts=seq(1,N,1) 
  Cp0=2/3 #Part of underflow model (see Follett et al. 2021)
  Cb=Cp0*coeff/S-1 #Part of underflow model (see Follett et al. 2021)
  
  for(i in Npts){ #For each channel section...
    colhj <- Hpts 
    
    #Jam present at section
    if (J[i]==1){
    colQNoJam <- B[i]*sqrt(9.8*colhj^3*S[i]/coeff[i]) 
    colQChannelBelowJam <- (colhj<a[i])*B[i]*sqrt(9.8*colhj^3*S[i]/coeff[i])
    
    #Case 1 (HJ<=HBF): weir over jam + FP resistance 
    if (Hj[i]<=HBF[i]){
      colQChannelJam <- (colhj<=Hj[i])*(colhj>=a[i])*B[i]*((2*9.8*(colhj-a[i])^3/(3^(3/2)*CA[i]))^(1/2)
                                           +(Cp0/(1+Cb[i]*a[i]/colhj)*9.8*a[i]^2*colhj)^(1/2))
      colQChannelJam[is.nan(colQChannelJam)] = 0
      
      colQChannelAboveJam <- (colhj>Hj[i])*(2/3*sqrt(2*9.8)*B[i]*(colhj-Hj[i])^(3/2)+max(colQChannelJam))
      matMin=matrix(,length(Hpts),2)
      matMin[,1]=colQNoJam
      matMin[,2]=colQChannelAboveJam
      
      colQChannelAboveJam <- (colhj>Hj[i])*rowMins(matMin,TRUE)
      #Q at an H can't be greater than uniform flow (drowned mode)
      colQChannelAboveJam[is.nan(colQChannelAboveJam)] = 0
      MaxQInChannel=max(c(max(colQChannelJam),max(colQChannelAboveJam*(colhj<=HBF[i]))))
      colQAboveBankfullFP <-(colhj>HBF[i])*((Beff[i]-1)*B[i])*((colhj-HBF[i])^(5/3)*sqrt(S[i])/nFP[i])
      colQAboveBankfullFP[is.nan(colQAboveBankfullFP)] = 0
      colQ=colQChannelBelowJam+colQChannelJam+colQChannelAboveJam+colQAboveBankfullFP
      
    #Case 2 (HJ>HBF): weir over jam and FP resistance  
    } else if (Hj[i]>HBF[i]){
#2. wide jam extending above bankfull depth (HJ>HBF)
    colQChannelJam <- (colhj<=Hj[i])*(colhj>=a[i])*B[i]*((2*9.8*(colhj-a[i])^3/(3^(3/2)*CA[i]))^(1/2)
                                                         +(Cp0/(1+Cb[i]*a[i]/colhj)*9.8*a[i]^2*colhj)^(1/2))
    colQChannelJam[is.nan(colQChannelJam)] = 0
  
    colQAboveBankfullJam <-(colhj>HBF[i])*(colhj<=Hj[i])*((Beff[i]-1)*B[i]*((2*9.8*(colhj-HBF[i])^3/(3^(3/2)*CA[i]))^(1/2)))
    colQAboveBankfullJam[is.nan(colQAboveBankfullJam)] = 0
    MaxQatHJ=max(colQChannelJam+colQAboveBankfullJam)

    colQAboveBankfullWeir <- (colhj>HBF[i])*(colhj>Hj[i])*(2*sqrt(2*9.8)/3*Beff[i]*B[i]*(colhj-Hj[i])^(3/2)+MaxQatHJ) #value when jam is full
    colQAboveBankfullWeir[is.nan(colQAboveBankfullWeir)] = 0
    
    colQ=colQChannelBelowJam+colQChannelJam+colQAboveBankfullJam+colQAboveBankfullWeir
    }
  #No jam
    } 
    else {
    colQChannel <- B[i]*sqrt(9.8*Hpts^3*S[i]/coeff[i])
    colQFP <- ((Beff[i]-1)*B[i])*((colhj-HBF[i])^(5/3)*sqrt(S[i])/nFP[i])
    colQFP[is.nan(colQFP)] = 0
    colQ<- colQChannel+colQFP
  }

    #find Q corresponding to H0 (unobstructed) for later interpolation
    colQH0Channel <- B[i]*sqrt(9.8*Hpts^3*S[i]/coeff[i])
    colQH0FP <-((Beff[i]-1)*B[i])*((colhj-HBF[i])^(5/3)*sqrt(S[i])/nFP[i])
    colQH0FP[is.nan(colQH0FP)] = 0
    
    colQH0<- colQH0Channel+colQH0FP
    colQH0 <- (colQH0>=colQ)*colQH0 + (colQH0<colQ)*colQ #drowned barrier
    colh0=interp1(colQH0,Hpts,colQ)
 
    colhj <- colhj*(colhj>colh0)+colh0*(colhj<=colh0)

    colVunif<-findunifvol(Ls[i],B[i],colh0,colhj,HBF[i],Beff[i])
    colVs<-colVunif
    
    if (J[i]==1){
      colVBWtotal<-findBWvolv2(Ls[i],B[i],colh0,colhj,HBF[i],Beff[i],S[i],colQ)

      colVs<-colVBWtotal+colVunif 
    }
    
    colV<-colVs
    colA <- colV/Ls[i] #longitudinal average cross-sectional area as in Hankin et al. 2020
    colLbw=(colhj-colh0)/S[i]
    rivarray[,,i] = c(colQ,colh0,colhj,colV,colA)

  }
  return(rivarray)
}
  
