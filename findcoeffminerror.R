
#choice of random or uniform jam distribution, site can be changed
#--------------
rivi=18; #index of site from field observations (default is Ouzel 2)
riv <- createriverRandomSpacing(0.65,allS[rivi],allB[rivi],allCA[rivi],allLs[rivi],allCf0[rivi],1)
#riv <- createriverUniformSpacing(0.65,allS[rivi],allB[rivi],allCA[rivi],allLs[rivi],allCf0[rivi],1)
#disp(c(allS[rivi],allB[rivi],allCA[rivi],allLs[rivi],allCf0[rivi]))
alloutJ=findoutflowv7(riv,coeff=1,jams=1)
tic()
#if changing site, need to change optimisation limits
coeffFit <- optimise(findcurverror,c(20,40),alloutJQ=alloutJ$Q,riv=riv,tol=0.1)
  toc()
#--------------------------
#calculates curves for unobstructed, to illustrate output
alloutN0<-findoutflowv7(riv,1,0)
plot(alloutJ$t_in/3600,alloutJ$Q_in,type='l',col='deepskyblue',xlab="time (h)",ylab=expression(paste("Discharge ","(m"^"3","/s)")),ylim=c(0.5*min(alloutJ$Q_in),max(1.2*alloutJ$Q_in)),xlim=c(0,1*72))
lines(alloutJ$t/3600,alloutJ$Q/1,col='deeppink3')
lines(alloutN0$t/3600,alloutN0$Q/1,col='orange1')
alloutNT<-findoutflowv7(riv,15,0)
lines(alloutNT$t/3600,alloutNT$Q/1,col='royalblue2')

alloutN<-findoutflowv7(riv,coeffFit$minimum,0)
lines(alloutN$t/3600,alloutN$Q/1,col='cadetblue')
alloutJR<-alloutJ

#calculate output parameters
tpeakOutJ=alloutJR$t[which.max(alloutJR$Q)]/3600
tpeakOutN0=alloutN0$t[which.max(alloutN0$Q)]/3600

QmaxIn=max(alloutJR$Q_in)
QmaxOutJ=max(alloutJR$Q)
QmaxOutN0=max(alloutN0$Q)

#display output
disp("Cfm QmaxIn,QmaxOutJ,QmaxOutN0,tpeakOutJ,tpeakOutN0")
disp(coeffFit$minimum,QmaxIn,QmaxOutJ,QmaxOutN0,tpeakOutJ,tpeakOutN0)

