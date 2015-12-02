# Plots for hypotheses using simply the means by species (original simulations)
## Created by RZF on Jul/27/15
library(grDevices)
col1<-rep(0,7)
col1[1]<-rgb(178,24,43, max = 255) # GD
col1[2]<-rgb(253,174,97, max = 255) # GST
col1[3]<-rgb(33,102,172, max = 255) #GST+GD
col1[4]<-rgb(239,138,98, max = 255) # MNT
col1[5]<-rgb(146,197,222,max=255)# MNT+GD
col1[6]<-rgb(171,171,171, max = 255)# GST+MNT
col1[7]<-rgb(0,0,0, max = 255) #GST+MNT+ GD
setwd('~/Dropbox/GenomeDownsizing/R programming/Simulations by Species')

#################GD only
bootstrapped.sim1<-read.table(file="bootstrap1.txt",sep="\n")[,1]
bootstrapped.sim<-matrix(bootstrapped.sim1,ncol=5, byrow=TRUE)
bootstrapped.sim<-bootstrapped.sim[-1,]
#Plot the ploidy vs. MSE estimate 
ret.rates<-seq(0.001,1,0.001) # retention rates tested
long1<-length(ret.rates)
total.colors<-rep(0,long1)
min.qdmean<-rep(0,long1) #Result of each hypothesis
for(j in 1:long1){
		min.qdmean[j]<-sum(bootstrapped.sim[j,])	
}

apply(bootstrapped.sim,2,which.min) #to obtain the best percentage at each ploidy under each hypothesis
apply(bootstrapped.sim,2,min) #best Total MSE
index.colors<-order(min.qdmean,decreasing=TRUE)
for(i in 1:5000){
  	total.colors[index.colors[i]]<-i/5000
  }
    col2<-rgb(0,0,0,total.colors)
    
plot(c(2,4,6,8,10,12),rep(0,6),type="l",ylim=c(0,100), col="white",xlab="Ploidy",ylab="Mean Squared Error Estimate")
for(j in 1:long1){
	points(c(4,6,8,10,12), bootstrapped.sim[j,], cex=3,type="l",col=col2[j])
}
besthyp1<-which.min(min.qdmean)
best.ratehyp1<-ret.rates[besthyp1] # 1- best.ratehyp 0.046
best.ateachploidy<-matrix(rep(0,5*7),nrow=7)
best.ateachploidy[1,]<-bootstrapped.sim[besthyp1,]
points(c(4,6,8,10,12), best.ateachploidy[1,], cex=5,type="l",col=col1[1],lwd=2)



############## GST and GST+GD
bootstrapped.sim3<-read.table(file="bootstrap5.txt",sep="\n")[,1]
bootstrapped.sim3<-matrix(bootstrapped.sim3,ncol=5, byrow=TRUE)
bootstrapped.sim3<-bootstrapped.sim3[-1,]
#Plot the ploidy vs. MSE estimate 
ret.rates<-seq(0.001,1,0.001) # retention rates tested
long1<-length(ret.rates)
total.colors<-rep(0,long1)
min.qdmean3<-rep(0,long1) #Result of each hypothesis
for(j in 1:long1){
		min.qdmean3[j]<-sum(bootstrapped.sim3[j,])	
}
apply(bootstrapped.sim3,2,which.min) #to obtain the best percentage at each ploidy under each hypothesis
apply(bootstrapped.sim3,2,min) #best Total MSE


index.colors<-order(min.qdmean3,decreasing=TRUE)
for(i in 1:5000){
  	total.colors[index.colors[i]]<-i/5000
  }
    col2<-rgb(0,0,0,total.colors)
    
plot(c(2,4,6,8,10,12),rep(0,6),type="l",ylim=c(0,100), col="white",xlab="Ploidy",ylab="Mean Squared Error Estimate")
for(j in 1:long1){
	points(c(4,6,8,10,12), bootstrapped.sim3[j,], cex=3,type="l",col=col2[j])
}
besthyp3<-which.min(min.qdmean3)
best.ratehyp3<-ret.rates[besthyp3] # 1- best.ratehyp 0.004
best.ateachploidy[3,]<-bootstrapped.sim3[besthyp3,]
points(c(4,6,8,10,12), best.ateachploidy[3,], cex=3,type="l",col=col1[3],lwd=2)

############## MNT and MNT+GD
bootstrapped.sim5<-read.table(file="mntgd.txt",sep="\n")[,1]
bootstrapped.sim5<-matrix(bootstrapped.sim5,ncol=5, byrow=TRUE)
bootstrapped.sim5<-bootstrapped.sim5[-1,]
#Plot the ploidy vs. MSE estimate 
ret.rates<-seq(0.001,1,0.001) # retention rates tested
long1<-length(ret.rates)
total.colors<-rep(0,long1)
min.qdmean5<-rep(0,long1) #Result of each hypothesis
for(j in 1:long1){
		min.qdmean5[j]<-sum(bootstrapped.sim5[j,])	
}
apply(bootstrapped.sim5,2,which.min) #to obtain the best percentage at each ploidy under each hypothesis
apply(bootstrapped.sim5,2,min) #best Total MSE



index.colors<-order(min.qdmean5,decreasing=TRUE)
for(i in 1:5000){
  	total.colors[index.colors[i]]<-i/5000
  }
    col2<-rgb(0,0,0,total.colors)
    
plot(c(2,4,6,8,10,12),rep(0,6),type="l",ylim=c(0,100), col="white",xlab="Ploidy",ylab="Mean Squared Error Estimate")
for(j in 1:long1){
	points(c(4,6,8,10,12), bootstrapped.sim5[j,], cex=3,type="l",col=col2[j])
}
besthyp5<-which.min(min.qdmean5)
best.ratehyp5<-ret.rates[besthyp5] # 1- best.ratehyp 0.058
best.ateachploidy[5,]<-bootstrapped.sim5[besthyp5,]
points(c(4,6,8,10,12), best.ateachploidy[5,], cex=3,type="l",col=col1[5],lwd=2)


############## GST+MNT and GST+MNT+GD
bootstrapped.sim7<-read.table(file="bootstrap6.txt",sep="\n")[,1]
bootstrapped.sim7<-matrix(bootstrapped.sim7,ncol=5, byrow=TRUE)
bootstrapped.sim7<-bootstrapped.sim7[-1,]
#Plot the ploidy vs. MSE estimate 
ret.rates<-seq(0.001,1,0.001) # retention rates tested
long1<-length(ret.rates)
total.colors<-rep(0,long1)
min.qdmean7<-rep(0,long1) #Result of each hypothesis
for(j in 1:long1){
		min.qdmean7[j]<-sum(bootstrapped.sim7[j,])	
}
apply(bootstrapped.sim7,2,which.min) #to obtain the best percentage at each ploidy under each hypothesis
apply(bootstrapped.sim7,2,min) #best Total MSE



index.colors<-order(min.qdmean7,decreasing=TRUE)
for(i in 1:5000){
  	total.colors[index.colors[i]]<-i/5000
  }
    col2<-rgb(0,0,0,total.colors)
    
plot(c(2,4,6,8,10,12),rep(0,6),type="l",ylim=c(0,100), col="white",xlab="Ploidy",ylab="Mean Squared Error Estimate")
for(j in 1:long1){
	points(c(4,6,8,10,12), bootstrapped.sim7[j,], cex=3,type="l",col=col2[j])
}
besthyp7<-which.min(min.qdmean7)
best.ratehyp7<-ret.rates[besthyp7] # 1- best.ratehyp 0.012
best.ateachploidy[7,]<-bootstrapped.sim7[besthyp7,]
points(c(4,6,8,10,12), best.ateachploidy[7,], cex=3,type="l",col=col1[7],lwd=2)

### All hypotheses
solutions<-read.table(file="bestateachploidy.txt",sep="\n")
solutions<-matrix(solutions$V1,ncol=5,byrow=TRUE)
best.ateachploidy[2,]<-solutions[3,]# GST
best.ateachploidy[4,]<-solutions[2,] #MNT
best.ateachploidy[6,]<-solutions[4,] #GST+MNT
par(mfrow=c(2,2))
plot(c(2,4,6,8,10,12),rep(0,6),type="l",ylim=c(0,35), col="white",xlab="Ploidy",ylab="Mean Squared Error Estimate")
for(j in 1:7){
	points(c(4,6,8,10,12), best.ateachploidy[j,], cex=3,type="l",col=col1[j],lwd=2)
}
legend(6,35,col=col1,lwd=2,legend=c("GD","GST","GST+GD","MNT","MNT+GD","GST+MNT","GST+MNT+GD"),cex=0.5)


###############
plot(1-ret.rates,min.qdmean,type="l",col=col1[1],xlab="Genome Downsizing Percentage",ylab="Total MSE",lwd=2)
points(1-best.ratehyp1,min(min.qdmean),pch=13, col=col1[1])

points(1-ret.rates,min.qdmean3,type="l",col=col1[3],lwd=2)
points(1-best.ratehyp3,min(min.qdmean3),pch=13, col=col1[3])
points(1-ret.rates,min.qdmean5,type="l",col=col1[5],lwd=2)
points(1-best.ratehyp5,min(min.qdmean5),pch=13, col=col1[5])
points(1-ret.rates,min.qdmean7,type="l",col=col1[7],lwd=2)
points(1-best.ratehyp7,min(min.qdmean7),pch=13, col=col1[7])

legend(0,280, col=col1[c(1,3,5,7)],lwd=2, legend=c("GD","GST+GD","MNT+GD","GST+MNT+GD"),cex=0.5)


######
min.GD<-min(min.qdmean)
rescale.GD<- exp(-min.qdmean+min.GD)
plot(1-ret.rates, rescale.GD,type="l" , xlab="Genome Downsizing Percentage", ylab="Approximate Profile Likelihood", col=col1[1])
min.GDGST<-min(min.qdmean3)
rescale.GDGST<-exp(-min.qdmean3+min.GDGST)
points(1-ret.rates, rescale.GDGST, type="l", col=col1[3])
min.GDMNT<-min(min.qdmean5)
rescale.GDMNT<-exp(-min.qdmean5+min.GDMNT)
points(1-ret.rates, rescale.GDMNT, type="l", col=col1[5])

min.GDGSTMNT<-min(min.qdmean7)
rescale.GDGSTMNT<-exp(-min.qdmean7+min.GDGSTMNT)
points(1-ret.rates, rescale.GDGSTMNT, type="l", col=col1[7])

legend(0.8,1, col=col1[c(1,3,5,7)],lwd=2, legend=c("GD","GST+GD","MNT+GD","GST+MNT+GD"),cex=0.6)

which(rescale.GD>.14042084)  # 651-652, 763-764  CI (0.236-0.349)
which(rescale.GDGST>.14042084) #773-774, 913-914 CI (0.086, 0.227)
which( rescale.GDMNT>.14042084)  #642-643 753-754  CI (0.247, 0.358)
which(rescale.GDGSTMNT>.14042084)#776-777, 917-918 CI (0.082, 0.224)

segments(0.236,.14042084, 0.358,.14042084,lty="dotted")
segments(0.236,-2, 0.236, .14042084,lty="dotted") 
segments(0.247,-2, 0.247, .14042084,lty="dotted") 
segments(0.349,-2, 0.349, .14042084,lty="dotted") 
segments(0.358,-2, 0.358, .14042084,lty="dotted") 

segments(0.082,.14042084, 0.227,.14042084,lty="dotted")
segments(0.082,-2, 0.082, .14042084,lty="dotted") 
segments(0.086,-2, 0.086, .14042084,lty="dotted") 
segments(0.227,-2, 0.227, .14042084,lty="dotted") 
segments(0.224,-2, 0.224, .14042084,lty="dotted") 


########################
########################
## Approximate AIC
AIC.GD=2+2*sum(best.ateachploidy[1,])
AIC.GST=2*5+2*sum(best.ateachploidy[2,])
AIC.GDGST=2*6+2*sum(best.ateachploidy[3,])
AIC.MNT=2*5+2*sum(best.ateachploidy[4,])
AIC.GDMNT=2*6+2*sum(best.ateachploidy[5,])
AIC.GSTMNT=2*5+2*sum(best.ateachploidy[6,])
AIC.GDGSTMNT=2*6+2*sum(best.ateachploidy[7,])