rm(list=ls())
# Polyploidy Series regression
# Rosana Zenil-Ferguson
# Last updated 

setwd('~/Dropbox/genomedownsizing/Datasets/')

#Reading Dataset for polyploidy series
angiocvalues=read.csv("ps2.csv")
#Defining parameters
complete<-which(angiocvalues$Chromploid==1)#268
total.ps<-which(angiocvalues$Minimum==1)
col2<-rgb(0,0,0,0.1)
#There are 
totalps<-max(angiocvalues[,13],na.rm=TRUE) #304

# 1-Family, 2-Genus, 3-Species, 4-Csome, 5-Ploidy, 6-Genome Size GS, 7-Repeated smaples, 8- Non multiple
#9- Minimum 10- Possible hybrid Mybehybrid 11- Both 12 Chromploid 13- PSID polyploidy series ID 
psploidylev<-angiocvalues$Ploidy
psgenome<-angiocvalues[,6]
#Relative differences
y.expected<-rep(0, 304)
y.observed<-rep(0,304)
obs.ploidy<-rep(0,304)
species.name<-rep(0,304)
# Plot differences in genome size
for(i in 1:304){
	aux<-which(angiocvalues$PSID==i)
	x.i<-psploidylev[aux]
	aux2<-x.i/x.i[1]
	y.i<-psgenome[aux]
	aux3<-angiocvalues$Name[aux]
	expected<-y.i[1]*aux2
	y.observed[i]<-y.i[2]
	y.expected[i]<- expected[2]
	obs.ploidy[i]<-x.i[2]
	species.name[i]<- as.character(aux3[1])
}	
aux4<-which(is.na(y.expected==TRUE))
y.expected<-y.expected[-aux4]
y.observed<-y.observed[-aux4]
obs.ploidy<-obs.ploidy[-aux4]
species.name<-species.name[-aux4]

table.autopolyploids<-data.frame(Name=species.name, GenomeObs=y.observed,GenomeExp=y.expected,Ploidy=obs.ploidy, row.names=species.name, stringsAsFactors=FALSE)

source("zanneangiospermtree.R")
library(MASS)
library(geiger)
library(ape)
library(nlme)
tmpploid<-treedata(angiosperm.tree, table.autopolyploids,sort=TRUE)
matched.datags<-tmpploid$data #121
matched.treegs<-tmpploid$phy
plot(matched.treegs,cex=0.3, type="fan")	

gs.obs<-as.numeric(matched.datags[,2])
gs.exp<-as.numeric(matched.datags[,3])
a<-cbind(gs.obs, gs.exp)
row.names(a)<-row.names(matched.datags)
write.table(a,"genomeregression.txt", sep=",",row.names=TRUE)
write.tree(matched.treegs,"autopolyploid.tre")
model.pagel1<-gls(gs.obs~gs.exp-1,correlation=corPagel(0.5,matched.treegs,fixed=FALSE), method="ML")

model.pagel0<-gls(gs.obs~gs.exp,correlation=corPagel(0,matched.treegs,fixed=TRUE), method="ML")

model.pagel0fixed<-gls(gs.obs~gs.exp-1,correlation=corPagel(0,matched.treegs,fixed=TRUE), method="ML")

model.OU1<-gls(gs.obs~gs.exp-1,correlation=corMartins(1,matched.treegs,fixed=FALSE), method="ML")

model.OU0<-gls(gs.obs~gs.exp-1,correlation=corMartins(0,matched.treegs,fixed=TRUE), method="ML")

model.BM<-gls(gs.obs~gs.exp-1,correlation=corBrownian(1,matched.treegs), method="ML")


##### Profile Likelihood for lambda
lambda0<-seq(0, 1, 0.005)
long1<-length(lambda0)
log.likes<-rep(0,long1)
for(i in 1:long1){
	model.i<-gls(gs.obs~gs.exp-1,correlation=corPagel(lambda0[i],matched.treegs,fixed=TRUE),method="ML")
	log.likes[i]<-logLik(model.i)[1]
}
rel.loglike<-exp(log.likes-(max(log.likes)))
plot(lambda0,rel.loglike, type="l", xlab=paste("Pagel's", expression(lambda),sep=" "),ylab="Relative profile likelihood")
##################
alpha0<-seq(0, 10, 0.005)
long1<-length(alpha0)
log.likes1<-rep(0,long1)
for(i in 1:long1){
	model.i<-gls(gs.obs~gs.exp-1,correlation=corMartins(alpha0[i],matched.treegs,fixed=TRUE),method="ML")
	log.likes1[i]<-logLik(model.i)[1]
}
rel.loglike1<-exp(log.likes1[-1]-(max(log.likes1[-1])))
plot(alpha0[-1],rel.loglike1, type="l", xlab=expression(alpha),ylab="Relative profile likelihood")

gamma0<-seq(0,10,0.005)

lk<-function(y,X,C, sig2){
    n<-nrow(C)
    v<-rep(0,n)
    V<-sig2*C+diag(v)
    beta<-solve(t(X)%*%solve(V)%*%X)%*%(t(X)%*%solve(V)%*%y)
    logL<--(1/2)*t(y-X%*%beta)%*%solve(V)%*%(y-X%*%beta)-
        (1/2)*determinant(V,logarithm=TRUE)$modulus[1]-
        (n/2)*log(2*pi)
     result<-cbind(beta,logL)
     return(result)
}

sigma0<-seq(0.05,2,0.005)
long1<-length(sigma0)
log.likes2<-rep(0, long1)
for(i in 1:long1){
	model.i<-lk(y=gs.obs, X=gs.exp, C=vcv(matched.treegs), sig2=sigma0[i])
	log.likes2[i]<-model.i[2]
	}
rel.loglike2<-exp(log.likes2-(max(log.likes2)))
plot(sigma0,rel.loglike2, type="l", xlab=expression(sigma^2),ylab="Relative profile likelihood",xlim=c(0,2))




$beta
[1] 0.7981314

$sig2e
[1] 1.000714

$logL
[1] -353.9077


6-2*(-353.9077)
713.8154















