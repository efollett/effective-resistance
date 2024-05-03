#takes input of inflow, channel network definition and returns outflow ("allout")
#calls subfunctions getinflow.R (defines inflow) 
#odes3.R (model system)

findoutflowv7 <- function(riv,coeff=1,jams=1){
  #coeff is multiple of channel unobstructed resistance!
  #jam=1 has jams as specified in createriver file; jams=0 is same channel with no jams
  tic()
library(pracma)  
library(deSolve)
library(Matrix)
  
#import river parameters from createriver() and inflow q,t from getinflow()
if (jams==1){
  #riv <- createriverExample()
  inflow <- getinflow(riv$minQ,riv$maxQ)
} else{
 # riv <- createriverExample()
 inflow <- getinflow(riv$minQ,riv$maxQ)
  riv$J <- 0*riv$J #same network as jams but with J=0 at all nodes
  riv$coeff <- coeff*riv$coeff #gives array with coeff=input in function call
}

#adds 1h to inflow time--avoids apparent bug where t queried by lsodes > max(inflow$t) 
inflow$t <- append(inflow$t,inflow$t[length(inflow$t)]+1)
inflow$Q <- append(inflow$Q,inflow$Q[length(inflow$Q)])

#set parameters for use in solving system
lx <- length(riv$x)+1 #used in later lines
N <- length(riv$x)-1; #number of segments--follows code accompanying Hankin et al. (2020)
#not sure why this isn't N=lx-1--perhaps to do with definition of inflow ex. as exiting first segment?
#ignores first and dummy segment

ad <- sparseMatrix(i=1:(N-1),j=2:N,x=1,dims=c(N,N)) #set up for 1 single channel; would need to alter for sub-channels in network
#for future improvement-we could compute ad, if the user inputs sub-channels & connection node...would need to specify order of solution 
l <- matrix(riv$x[2:(lx-1)]-riv$x[1:(lx-2)],N,1)
S <- matrix(-(riv$z[2:(lx-1)]-riv$z[1:(lx-2)])/(riv$x[2:(lx-1)]-riv$x[1:(lx-2)]),N,1)
CA <- matrix(riv$CA[2:(lx-1)],N,1) 
g <- matrix(rep(9.8,N),N,1)
B <- matrix(riv$B[2:(lx-1)],N,1)
coeff <- matrix(riv$coeff[2:(lx-1)],N,1)
J <- matrix(riv$J[2:(lx-1)],N,1)
a <- matrix(riv$a[2:(lx-1)],N,1)
HJ <- matrix(riv$HJ[2:(lx-1)],N,1)
HBF <- matrix(riv$HBF[2:(lx-1)],N,1)
Beff <- matrix(riv$Beff[2:(lx-1)],N,1)
nFP <- matrix(riv$nFP[2:(lx-1)],N,1)

#Calculate precalculated hydraulic table 
rivarray <-makelookuparrayv6(CA,S,coeff,B,l,N,J,a,HJ,HBF,Beff,riv$QBF,nFP)
              
#solve system
#1. constant inflow: find initial conditions (ICs)
parmsIC <- list(B,S,coeff,g,CA,J,inflow$t,rep(inflow$Q[1],length(inflow$t)),N,ad,l,a,HJ,rivarray)
baseflowguess=(B[1]*(riv$coeff[1]*(min(inflow$Q)/B[1])^2/(g[1]*S[1]))^(1/3)) #guess area
yIC=rep(baseflowguess,N) 
#disp(baseflowguess)
outIC <- lsodes(y=yIC,times=c(0.25,3*72*3600),func=odes3,parms=parmsIC)
#2. variable inflow: use ICs from (1)
parms<-list(B,S,coeff,g,CA,J,inflow$t,inflow$Q,N,ad,l,a,HJ,rivarray)
times=seq(0, (3*72)*3600,0.25*3600) #related to bug mentioned in Line 12-input only original time which reduces maximum t input by lsodes
outA <- lsodes(y=as.vector(outIC[2,2:(lx-1)]),times=times,func=odes3,parms=parms)
outA <- outA[,2:(lx-1)] #in R first entry is time

#find h,Q only for last node, N, to reduce time
#save in allout structure
allout <- list(A = outA[,N])
lenTimes=length(times)
hOut=rep(0,lenTimes)
h0Out=rep(0,lenTimes)
hJOut=rep(0,lenTimes)
QOut=rep(0,lenTimes)

for(i in seq(1,lenTimes)){
  colchoice<-J[N]+2
  hOut[i]=interp1(rivarray[,5,N],rivarray[,colchoice,N],outA[i,N])
  QOut[i]=interp1(rivarray[,5,N],rivarray[,1,N],outA[i,N])
  h0Out[i]=interp1(rivarray[,5,N],rivarray[,2,N],outA[i,N])
  hJOut[i]=interp1(rivarray[,5,N],rivarray[,3,N],outA[i,N])
}

allout["h"] <- list(h = hOut)
allout["h0"] <- list(h = h0Out)
allout["hJ"] <- list(h = hJOut)
allout["Q"] <- list(Q = QOut)
allout["t"] <- list(t = times)
allout["Q_in"] <- list(Q_in = inflow$Q)
allout["t_in"] <- list(t_in = inflow$t*3600)

disp('model run complete')
toc() #paired with tic(); finds how long this took to run
return(allout)
}
