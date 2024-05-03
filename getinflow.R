#computes inflow curve for model
getinflow <- function(minQ,maxQ){
mu = 6;
sigma = 1;

#Compute wave
 xi=seq(0,72*10,0.05) #hr
yi=1/(sigma*sqrt(2*pi))*exp(-0.5*((xi-mu)/sigma)^2)
coeff=(maxQ-minQ)/max(yi);
yi=yi*coeff+minQ;

#create list for output
inflow <- list(t=xi) #h-for ease of entry, but could change to s for correspondence with Q 
inflow["Q"] <- list(Q=yi) #m3/s
return(inflow)
}

