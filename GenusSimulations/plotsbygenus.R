library(grDevices)
col1<-rep(0,7)
col1[1]<-rgb(178,24,43, max = 255) # GD
col1[2]<-rgb(253,174,97, max = 255) # GST
col1[3]<-rgb(33,102,172, max = 255) #GST+GD
col1[4]<-rgb(239,138,98, max = 255) # MNT
col1[5]<-rgb(146,197,222,max=255)# MNT+GD
col1[6]<-rgb(171,171,171, max = 255)# GST+MNT
col1[7]<-rgb(0,0,0, max = 255) #GST+MNT+ GD
setwd("~/Dropbox/GenomeDownsizing/R programming/Simulations by Genus/")

#################GD only
bootstrapped.sim1<-read.table(file="bootstrapbygenusgdnewaverage.txt",sep="\n")[,1]
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

# apply(bootstrapped.sim,2,which.min) to obtain the best percentage at each ploidy under each hypothesis
#apply(bootstrapped.sim,2,min) best Total MSE

[1]  0.006508806  5.125497581  0.937028077 19.560964868  0.713035180
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
best.ratehyp1<-ret.rates[besthyp1] # 1- best.ratehyp 0.003
best.ateachploidy<-matrix(rep(0,5*7),nrow=7)
best.ateachploidy[1,]<-bootstrapped.sim[besthyp1,]
min(min.qdmean) #29.51766
points(c(4,6,8,10,12), best.ateachploidy[1,], cex=5,type="l",col=col1[1],lwd=2)

############## GST and GST+GD
bootstrapped.sim2<-read.table(file="bootstrapbygenusgs.txt",sep="\n")[,1]
#points(c(4,6,8,10,12), bootstrapped.sim2, cex=3,type="l",col=col1[2],lwd=2)
best.ateachploidy[2,]<-bootstrapped.sim2
sum(bootstrapped.sim2)

bootstrapped.sim3<-read.table(file="bootstrapbygenusgsgd.txt",sep="\n")[,1]
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
#apply(bootstrapped.sim3,2,which.min) # to obtain the best percentage at each ploidy under each hypothesis
apply(bootstrapped.sim,2,min) #best Total MSE


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
min(min.qdmean3)
best.ratehyp3<-ret.rates[besthyp3] # 1- best.ratehyp 0.000
best.ateachploidy[3,]<-bootstrapped.sim3[besthyp3,]

points(c(4,6,8,10,12), best.ateachploidy[3,], cex=5,type="l",col=col1[3],lwd=2)

############## MNT and MNT+GD
bootstrapped.sim4<-read.table(file="bootstrapbygenusmn.txt",sep="\n")[,1]
plot(c(4,6,8,10,12), bootstrapped.sim4, cex=3,type="l",col=col1[4],lwd=2)
best.ateachploidy[4,]<-bootstrapped.sim4

bootstrapped.sim5<-read.table(file="bootstrapbygenusmngd.txt",sep="\n")[,1]
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
apply(bootstrapped.sim5,2,which.min) # to obtain the best percentage at each ploidy under each hypothesis
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
min(min.qdmean5)
best.ratehyp5<-ret.rates[besthyp5] # 1- best.ratehyp 0.058
best.ateachploidy[5,]<-bootstrapped.sim5[besthyp5,]
points(c(4,6,8,10,12), best.ateachploidy[5,], cex=5,type="l",col=col1[5],lwd=2)


############## GST+MNT and GST+MNT+GD
bootstrapped.sim6<-read.table(file="bootstrapbygenusgsmn.txt",sep="\n")[,1]
#plot(c(4,6,8,10,12), bootstrapped.sim6, cex=3,type="l",col=col1[6],lwd=2)
best.ateachploidy[6,]<-bootstrapped.sim6

bootstrapped.sim7<-read.table(file="bootstrapbygenusgsmngd.txt",sep="\n")[,1]
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
apply(bootstrapped.sim7,2,which.min) # to obtain the best percentage at each ploidy under each hypothesis
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
sum(best.ateachploidy[7,])
points(c(4,6,8,10,12), best.ateachploidy[7,], cex=5,type="l",col=col1[7],lwd=2)

### All hypotheses
plot(c(2,4,6,8,10,12),rep(0,6),type="l",ylim=c(0,30), col="white",xlab="Ploidy",ylab="Mean Squared Error Estimate")
for(j in 1:7){
	points(c(4,6,8,10,12), best.ateachploidy[j,], cex=3,type="l",col=col1[j],lwd=2)
}
legend(3,25,col=col1,lwd=2,legend=c("GD","GST","GST+GD","MNT","MNT+GD","GST+MNT","GST+MNT+GD"))


###############
plot(1-ret.rates,min.qdmean,type="l",col=col1[1],xlab="Genome Downsizing Percentage",ylab="Total MSE",lwd=2)
points(1-best.ratehyp1,min(min.qdmean),pch=13, col=col1[1])

points(1-ret.rates,min.qdmean3,type="l",col=col1[3],lwd=2)
points(1-best.ratehyp3,min(min.qdmean3),pch=13, col=col1[3])
points(1-ret.rates,min.qdmean5,type="l",col=col1[5],lwd=2)
points(1-best.ratehyp5,min(min.qdmean5),pch=13, col=col1[5])
points(1-ret.rates,min.qdmean7,type="l",col=col1[7],lwd=2)
points(1-best.ratehyp7,min(min.qdmean7),pch=13, col=col1[7])

legend(0.2,300, col=col1[c(1,3,5,7)],lwd=2, legend=c("GD","GD+GST","GD+MNT","GD+GST+MNT"))

#####
min.GD<-min(min.qdmean)
rescale.GD<- exp(-min.qdmean+min.GD)
plot(1-ret.rates, rescale.GD,type="l" , xlab="Genome Downsizing Percentage", ylab="Approximate Profile Likelihood", col=col1[1],xlim=c(0,0.14),lwd=2)
min.GDGST<-min(min.qdmean3)
rescale.GDGST<-exp(-min.qdmean3+min.GDGST)
points(1-ret.rates, rescale.GDGST, type="l", col=col1[3],lwd=2)
min.GDMNT<-min(min.qdmean5)
rescale.GDMNT<-exp(-min.qdmean5+min.GDMNT)
points(1-ret.rates, rescale.GDMNT, type="l", col=col1[5],lwd=2)

min.GDGSTMNT<-min(min.qdmean7)
rescale.GDGSTMNT<-exp(-min.qdmean7+min.GDGSTMNT)
points(1-ret.rates, rescale.GDGSTMNT, type="l", col=col1[7],lwd=2)

legend(0.1,1, col=col1[c(1,3,5,7)],lwd=2, legend=c("GD","GST+GD","MNT+GD","GST+MNT+GD"),cex=0.6)

which(rescale.GD>.14042084)  # 957-958,   100  CI (0-0.043)
which(rescale.GDGST>.14042084) #, 984-985, 100 CI (0, 0.015)
which( rescale.GDMNT>.14042084)  #952-953 1000  CI (0, 0.048)
which(rescale.GDGSTMNT>.14042084)#982-983, 100 CI (0, 0.018)

segments(0,.14042084, 0.048,.14042084,lty="dotted")
segments(0,-2, 0, .14042084,lty="dotted") 
segments(0.043,-2, 0.043, .14042084,lty="dotted") 
segments(0.015,-2, 0.015, .14042084,lty="dotted") 
segments(0.048,-2, 0.048, .14042084,lty="dotted") 
segments(0.018,-2, 0.018, .14042084,lty="dotted") 

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

