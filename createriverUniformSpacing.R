#creates 1D river channel with uniformly distributed jams
createriverUniformSpacing <-function(HBF,S,BBF,CA,Ls,Cf0,Beffi){
g=9.8
#HBF=0.65 #m const. for all reaches

#calculate inflow
QBF=BBF*sqrt(g*HBF^3*S/Cf0) #check if close, then use this one
QBFJ=BBF*sqrt(2*g*HBF^3/(CA*3^(3/2)))
maxQ=QBFJ
minQ=0.5*QBFJ

N=100; #Number of segments

LR = Ls*(N+1); #because Nsegs=Lx-2 in model
Nsegs = N
slope=S
Ns=N+2 #first segment, plus dummy segment at end for uniform flow

riv <- list(CA=c(0,rep(CA,N),0))
riv["J"] <-list(J=c(0,rep(1,N),0))
riv["x"] <- list(x=c(0,seq(Ls,(N*Ls),Ls))) 
riv["z"] <- list(z=S*max(riv$x)-S*riv$x+0)
riv["B"] <- list(B=rep(BBF,length(riv$x))) #m section bankfull width 
riv["S"] <- list(slope=rep(S,length(riv$x))) #[] section slope
riv["QBF"] <- QBF #m3/s average bankfull discharge for river section--used by getinflow.R--only one number as inflow used for whole reach
riv["coeff"] <- list(coeff=rep(Cf0,length(riv$x))) #[] unobstructed bed resistance
riv["a"] <- list(a=rep(0,length(riv$x))) #m vertical width of lower gap 
riv["HJ"] <- list(HJ=rep(100*HBF,length(riv$x))) #distance of top edge of jam above bed
riv["HBF"] <-list(HBF=rep(HBF,length(riv$x))) 
riv["Beff"]<-list(Beff=rep(Beffi,length(riv$x))) #Beff*BBF=channel width above bankfull depth; Beff=1 has no "floodplain"
riv["nFP"]<-list(nFP=rep(0.1,length(riv$x))) #manning's floodplain resistance (due to Yochum document for estimation)
riv["minQ"]<-minQ
riv["maxQ"]<-maxQ
riv["LR"]<-LR

return(riv)
}

