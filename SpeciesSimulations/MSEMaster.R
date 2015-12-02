# Bootstrap simulations for second type of selection
# Rosana Zenil-Ferguson
#Last updated 11/18/14

# Three libraries to manipulate trees
setwd("/Users/Roszenil/Dropbox/Summer 2013/DescriptiveStats")
library(grDevices)
col1<-rep(0,6)
col1[1]<-rgb(136,204,238, max = 255)
col1[2]<-rgb(17,119,51, max = 255)
col1[3]<-rgb(221,204,119, max = 255)
col1[4]<-rgb(204,102,119, max = 255)
col1[5]<-rgb(170,68,153, max = 255)
col1[6]<-rgb(51,34,136, max = 255)

#Reading Dataset 
angiocvalues=read.csv("apgIIcurated2.csv")
#Defining parameters
complete<-which(angiocvalues$Chromploid==1)#5030
just.chrom<-which(angiocvalues$Minimum==1)#7124
polyploidseries<-angiocvalues[,15]

# 1-Family, 2-Genus, 3-Species, 4-Authority, 5-Chromosome.number, 6-Ploidy.level, 7-Estimation.method, 8-Voucher, 9-X1C..pg., 10-Original.Reference,
# 11-Lifecycle, 12-Group, 13-Paper, 14-Name, 15-History.of.Polyploidy, 16-Repeated.Sample, 17-Non.multiple, 18-Minimum, 19-Missing.Data, 20-Possible.Hybrid, 
#21-Both, 22-Chromploid
 

#which.genus<-as.vector(unlist(levels(genus)))
#how.many.genus<-length(which.genus)

n.complete<-angiocvalues[complete,5]/angiocvalues[complete,6]
fam.complete<-angiocvalues[complete,1]#248 (all of them)
genus.complete<-angiocvalues[complete,2]
species.complete<-angiocvalues[complete,3]
ploidy.complete<-angiocvalues[complete,6]

#######

genome.size<-angiocvalues[just.chrom,9]
ploidy.genome<-angiocvalues[just.chrom,6]
mono.number<-angiocvalues[just.chrom,5]/angiocvalues[just.chrom,6]
genome.boxplot<-boxplot(genome.size~ploidy.genome)



### How many of each removing repeated samples and polyploidy series considering only the lowest ploidy
diploids<-which(ploidy.genome==2)
triploids<-which(ploidy.genome==3)
tetraploids<-which(ploidy.genome==4)
pentaploids<-which(ploidy.genome==5)
hexaploids<-which(ploidy.genome==6)
heptaploids<-which(ploidy.genome==7)
octoploids<-which(ploidy.genome==8)
enneaploids<-which(ploidy.genome==9)
decaploids<-which(ploidy.genome==10)
undecaploids<-which(ploidy.genome==11)
duodecaploids<-which(ploidy.genome==12)

how.many2x<-length(diploids)
how.many3x<-length(triploids)
how.many4x<-length(tetraploids)
how.many5x<-length(pentaploids)
how.many6x<-length(hexaploids)
how.many7x<-length(heptaploids)
how.many8x<-length(octoploids)
how.many9x<-length(enneaploids)
how.many10x<-length(decaploids)
how.many11x<-length(undecaploids)
how.many12x<-length(duodecaploids)



#### Empirical quantiles of the threshold
ecdf(genome.size[diploids])(53.4)
ecdf(genome.size[triploids])(77.55)
ecdf(genome.size[tetraploids])(80.6)
ecdf(genome.size[pentaploids])(26.65)
ecdf(genome.size[hexaploids])(78.7)
ecdf(genome.size[heptaploids])(6.93)
ecdf(genome.size[octoploids])(47.48)
ecdf(genome.size[enneaploids])(20.28)
ecdf(genome.size[decaploids])(12.29)
ecdf(genome.size[undecaploids])(1.12)
ecdf(genome.size[duodecaploids])(12.84)

####HYPOTHESIS 1 JUST GENOME DOWNSIZING
# Function where retention rates are given rr and how many bootstrap samples
retention.rate<-function(rr,B){ 
	#Bootstrapp means and standard error
bootstrapped.mean3x<-rep(0,B)
bootstrapped.sd3x<-rep(0,B)
bootstrapped.mean4x<-rep(0,B)
bootstrapped.sd4x<-rep(0,B)
bootstrapped.mean5x<-rep(0,B)
bootstrapped.sd5x<-rep(0,B)
bootstrapped.mean6x<-rep(0,B)
bootstrapped.sd6x<-rep(0,B)
bootstrapped.mean7x<-rep(0,B)
bootstrapped.sd7x<-rep(0,B)
bootstrapped.mean8x<-rep(0,B)
bootstrapped.sd8x<-rep(0,B)
bootstrapped.mean9x<-rep(0,B)
bootstrapped.sd9x<-rep(0,B)
bootstrapped.mean10x<-rep(0,B)
bootstrapped.sd10x<-rep(0,B)
bootstrapped.mean11x<-rep(0,B)
bootstrapped.sd11x<-rep(0,B)
bootstrapped.mean12x<-rep(0,B)
bootstrapped.sd12x<-rep(0,B)


 #retention rate
for (i in 1:B){
	
	genome.sample<-genome.size[diploids]
		
	#Tetraploids
bootstrap.sample4x<-sample(x=genome.sample, size=how.many4x, replace = TRUE, prob = NULL)*2*rr
bootstrapped.mean4x[i]<-mean(bootstrap.sample4x) # Mean of the new bootstrap sample
bootstrapped.sd4x[i]<-sd(bootstrap.sample4x) # Standard error of the new bootstrap sample
	
	
	genome.sample<-genome.size[tetraploids]
	
	#Hexaploids
bootstrap.sample6x<-sample(x=genome.sample, size=how.many6x, replace = TRUE, prob = NULL)*1.5*rr
bootstrapped.mean6x[i]<-mean(bootstrap.sample6x)
bootstrapped.sd6x[i]<-sd(bootstrap.sample6x)


	genome.sample<-genome.size[hexaploids]
#Octoploids
bootstrap.sample8x<-sample(x=genome.sample,size=how.many8x,replace=TRUE,prob=NULL)*1.33*rr
bootstrapped.mean8x[i]<-mean(bootstrap.sample8x)
bootstrapped.sd8x[i]<-sd(bootstrap.sample8x)


	genome.sample<-genome.size[octoploids]
	
	#Decaploids
bootstrap.sample10x<-sample(x=genome.sample, size=how.many10x, replace = TRUE, prob = NULL)*1.25*rr
bootstrapped.mean10x[i]<-mean(bootstrap.sample10x)
bootstrapped.sd10x[i]<-sd(bootstrap.sample10x)

 
	genome.sample<-genome.size[decaploids]
	

bootstrap.sample12x<-sample(x=genome.sample, size=how.many12x, replace = TRUE, prob = NULL)*1.2*rr
bootstrapped.mean12x[i]<-mean(bootstrap.sample12x)
bootstrapped.sd12x[i]<-sd(bootstrap.sample12x)
}

genome.means<-c(mean(genome.size[diploids]),mean(genome.size[tetraploids]),mean(genome.size[hexaploids]) ,mean(genome.size[octoploids]),mean(genome.size[decaploids]),mean(genome.size[duodecaploids])) # observed means for each ploidy


genome.sd<-c(sd(genome.size[diploids]),sd(genome.size[tetraploids]),sd(genome.size[hexaploids]),sd(genome.size[octoploids]),sd(genome.size[decaploids]),sd(genome.size[duodecaploids])) #observed standard error for each ploidy

# Mean Squared error- Distance from observed vs bootstraped the closer to 0 the better is the prediction
qdmean4x<-sum((bootstrapped.mean4x-genome.means[2])^2)/B
qdmean6x<-sum((bootstrapped.mean6x-genome.means[3])^2)/B
qdmean8x<-sum((bootstrapped.mean8x-genome.means[4])^2)/B
qdmean10x<-sum((bootstrapped.mean10x-genome.means[5])^2)/B
qdmean12x<-sum((bootstrapped.mean12x-genome.means[6])^2)/B

qdmeanh1<-c(qdmean4x,qdmean6x,qdmean8x,qdmean10x,qdmean12x) # Vector of distances
return(qdmeanh1)
}


B<-10000 #number of bootstraps
ret.rates<-seq(0.001,1,0.001) # retention rates tested
long1<-length(ret.rates)
min.qdmean<-rep(0,long1) #Result of each hypothesis
bootstrapped.sim1<-rep(0,5)# Samples with each rate
for(j in 1:long1){
	qdmean<-retention.rate(rr=ret.rates[j],B)
	bootstrapped.sim1<-rbind(bootstrapped.sim1,qdmean)
	min.qdmean[j]<-sum(qdmean)	
}

write.table(bootstrapped.sim1,file="bootstrap1.txt",sep="\n",row.names=FALSE,col.names=FALSE)
#Plot the ploidy vs. MSE estimate 
total.colors<-rep(0,long1)
index.colors<-order(min.qdmean,decreasing=TRUE)
for(i in 1:5000){
  	total.colors[index.colors[i]]<-i/5000
  }
    col2<-rgb(0,0,0,total.colors)


plot(c(2,4,6,8,10,12),rep(0,6),type="l",ylim=c(0,100), col="white",xlab="Ploidy",ylab="Mean Squared Error Estimate")
for(j in 2:long1){
	points(c(4,6,8,10,12), bootstrapped.sim1[j,], cex=3,type="l",col=col2[j])
}
besthyp1<-which.min(min.qdmean)
best.ratehyp1<-ret.rates[besthyp1] # 1- best.ratehyp 0.301
best.ateachploidy<-matrix(rep(0,5*6),nrow=6)
best.ateachploidy[1,]<-bootstrapped.sim1[besthyp1+1,]
points(c(4,6,8,10,12), best.ateachploidy[1,], cex=3,type="l",col=col1[1],lwd=2)



#####################
##############################################################################################################
##############HYPOTHESIS 2: JUST MONOPLOID NUMBER THRESHOLD
retention.rate2<-function(rr,B){

bootstrapped.mean3x<-rep(0,B)
bootstrapped.sd3x<-rep(0,B)
bootstrapped.mean4x<-rep(0,B)
bootstrapped.sd4x<-rep(0,B)
bootstrapped.mean5x<-rep(0,B)
bootstrapped.sd5x<-rep(0,B)
bootstrapped.mean6x<-rep(0,B)
bootstrapped.sd6x<-rep(0,B)
bootstrapped.mean7x<-rep(0,B)
bootstrapped.sd7x<-rep(0,B)
bootstrapped.mean8x<-rep(0,B)
bootstrapped.sd8x<-rep(0,B)
bootstrapped.mean9x<-rep(0,B)
bootstrapped.sd9x<-rep(0,B)
bootstrapped.mean10x<-rep(0,B)
bootstrapped.sd10x<-rep(0,B)
bootstrapped.mean11x<-rep(0,B)
bootstrapped.sd11x<-rep(0,B)
bootstrapped.mean12x<-rep(0,B)
bootstrapped.sd12x<-rep(0,B)



#retention rate
for (i in 1:B){
	
	genome.sample<-genome.size[diploids]
	mono.sample<-mono.number[diploids]
	individualstoduplicate<-which(mono.sample<31)	
		
	#Tetraploids
bootstrap.sample4x<-sample(x=genome.sample[individualstoduplicate], size=how.many4x, replace = TRUE, prob = NULL)*2*rr
bootstrapped.mean4x[i]<-mean(bootstrap.sample4x)
bootstrapped.sd4x[i]<-sd(bootstrap.sample4x)
	
	
	genome.sample<-genome.size[tetraploids]
	mono.sample<-mono.number[tetraploids]
	individualstoduplicate<-which(mono.sample<31)	
	
	#Hexaploids
bootstrap.sample6x<-sample(x=genome.sample[individualstoduplicate], size=how.many6x, replace = TRUE, prob = NULL)*1.5*rr
bootstrapped.mean6x[i]<-mean(bootstrap.sample6x)
bootstrapped.sd6x[i]<-sd(bootstrap.sample6x)


	genome.sample<-genome.size[hexaploids]
	mono.sample<-mono.number[hexaploids]
	individualstoduplicate<-which(mono.sample<31)
#Octoploids
bootstrap.sample8x<-sample(x=genome.sample[individualstoduplicate],size=how.many8x,replace=TRUE,prob=NULL)*1.33*rr
bootstrapped.mean8x[i]<-mean(bootstrap.sample8x)
bootstrapped.sd8x[i]<-sd(bootstrap.sample8x)


	genome.sample<-genome.size[octoploids]
	mono.sample<-mono.number[octoploids]
	individualstoduplicate<-which(mono.sample<25)
	
	#Decaploids
bootstrap.sample10x<-sample(x=genome.sample[individualstoduplicate], size=how.many10x, replace = TRUE, prob = NULL)*1.25*rr
bootstrapped.mean10x[i]<-mean(bootstrap.sample10x)
bootstrapped.sd10x[i]<-sd(bootstrap.sample10x)

 
	genome.sample<-genome.size[decaploids]
	mono.sample<-mono.number[decaploids]
	individualstoduplicate<-which(mono.sample<10)
	

bootstrap.sample12x<-sample(x=genome.sample[individualstoduplicate], size=how.many12x, replace = TRUE, prob = NULL)*1.2*rr
bootstrapped.mean12x[i]<-mean(bootstrap.sample12x)
bootstrapped.sd12x[i]<-sd(bootstrap.sample12x)
}

genome.means<-c(mean(genome.size[diploids]),mean(genome.size[tetraploids]),mean(genome.size[hexaploids]) ,mean(genome.size[octoploids]),mean(genome.size[decaploids]),mean(genome.size[duodecaploids]))

genome.sd<-c(sd(genome.size[diploids]),sd(genome.size[tetraploids]),sd(genome.size[hexaploids]),sd(genome.size[octoploids]),sd(genome.size[decaploids]),sd(genome.size[duodecaploids]))

qdmean4x<-sum((bootstrapped.mean4x-genome.means[2])^2)/B
qdmean6x<-sum((bootstrapped.mean6x-genome.means[3])^2)/B
qdmean8x<-sum((bootstrapped.mean8x-genome.means[4])^2)/B
qdmean10x<-sum((bootstrapped.mean10x-genome.means[5])^2)/B
qdmean12x<-sum((bootstrapped.mean12x-genome.means[6])^2)/B


qdmeanh1<-c(qdmean4x,qdmean6x,qdmean8x,qdmean10x,qdmean12x)
return(qdmeanh1)
}



qdmean<-retention.rate2(rr=1,B)
best.ateachploidy[2,]<-qdmean



points(c(4,6,8,10,12), best.ateachploidy[2,], cex=3,type="l",col=col1[2],lwd=2)




#####################
##############################################################################################################
##############HYPOTHESIS 3: JUST GENOME SIZE THRESHOLD
retention.rate3<-function(rr,B){

bootstrapped.mean3x<-rep(0,B)
bootstrapped.sd3x<-rep(0,B)
bootstrapped.mean4x<-rep(0,B)
bootstrapped.sd4x<-rep(0,B)
bootstrapped.mean5x<-rep(0,B)
bootstrapped.sd5x<-rep(0,B)
bootstrapped.mean6x<-rep(0,B)
bootstrapped.sd6x<-rep(0,B)
bootstrapped.mean7x<-rep(0,B)
bootstrapped.sd7x<-rep(0,B)
bootstrapped.mean8x<-rep(0,B)
bootstrapped.sd8x<-rep(0,B)
bootstrapped.mean9x<-rep(0,B)
bootstrapped.sd9x<-rep(0,B)
bootstrapped.mean10x<-rep(0,B)
bootstrapped.sd10x<-rep(0,B)
bootstrapped.mean11x<-rep(0,B)
bootstrapped.sd11x<-rep(0,B)
bootstrapped.mean12x<-rep(0,B)
bootstrapped.sd12x<-rep(0,B)

#retention rate
for (i in 1:B){
	
	genome.sample<-genome.size[diploids]
	mono.sample<-mono.number[diploids]
	individualstoduplicate<-which(genome.sample< 53.4)	
		
	#Tetraploids
bootstrap.sample4x<-sample(x=genome.sample[individualstoduplicate], size=how.many4x, replace = TRUE, prob = NULL)*2*rr
bootstrapped.mean4x[i]<-mean(bootstrap.sample4x)
bootstrapped.sd4x[i]<-sd(bootstrap.sample4x)
	
	
	genome.sample<-genome.size[tetraploids]
	mono.sample<-mono.number[tetraploids]
	individualstoduplicate<-which(genome.sample<80.61)	
	
	#Hexaploids
bootstrap.sample6x<-sample(x=genome.sample[individualstoduplicate], size=how.many6x, replace = TRUE, prob = NULL)*1.5*rr
bootstrapped.mean6x[i]<-mean(bootstrap.sample6x)
bootstrapped.sd6x[i]<-sd(bootstrap.sample6x)


	genome.sample<-genome.size[hexaploids]
	mono.sample<-mono.number[hexaploids]
	individualstoduplicate<-which(genome.sample<78.71)
#Octoploids
bootstrap.sample8x<-sample(x=genome.sample[individualstoduplicate],size=how.many8x,replace=TRUE,prob=NULL)*1.33*rr
bootstrapped.mean8x[i]<-mean(bootstrap.sample8x)
bootstrapped.sd8x[i]<-sd(bootstrap.sample8x)


	genome.sample<-genome.size[octoploids]
	mono.sample<-mono.number[octoploids]
	individualstoduplicate<-which(genome.sample< 47.49)
	
	#Decaploids
bootstrap.sample10x<-sample(x=genome.sample[individualstoduplicate], size=how.many10x, replace = TRUE, prob = NULL)*1.25*rr
bootstrapped.mean10x[i]<-mean(bootstrap.sample10x)
bootstrapped.sd10x[i]<-sd(bootstrap.sample10x)

 
	genome.sample<-genome.size[decaploids]
	mono.sample<-mono.number[decaploids]
	individualstoduplicate<-which(genome.sample< 12.30)
	

bootstrap.sample12x<-sample(x=genome.sample[individualstoduplicate], size=how.many12x, replace = TRUE, prob = NULL)*1.2*rr
bootstrapped.mean12x[i]<-mean(bootstrap.sample12x)
bootstrapped.sd12x[i]<-sd(bootstrap.sample12x)
}

genome.means<-c(mean(genome.size[diploids]),mean(genome.size[tetraploids]),mean(genome.size[hexaploids]) ,mean(genome.size[octoploids]),mean(genome.size[decaploids]),mean(genome.size[duodecaploids]))

# plot(c(2,4,6,8,10,12),genome.means,xlim=c(1,14),pch=19,ylim=c(0,30))
# points(c(2,4,6,8,10,12),genome.means,type='l')
# points(rep(4,100),bootstrapped.mean4x,col="red")
# points(rep(6,100),bootstrapped.mean6x,col="red")
# points(rep(8,100),bootstrapped.mean8x,col="red")
# points(rep(10,100),bootstrapped.mean10x,col="red")
# points(rep(12,100),bootstrapped.mean12x,col="red")


genome.sd<-c(sd(genome.size[diploids]),sd(genome.size[tetraploids]),sd(genome.size[hexaploids]),sd(genome.size[octoploids]),sd(genome.size[decaploids]),sd(genome.size[duodecaploids]))

# plot(c(2,4,6,8,10,12),genome.sd,xlim=c(1,14),pch=19,ylim=c(0,30))
# points(c(2,4,6,8,10,12),genome.sd,type='l')
# points(rep(4,100),bootstrapped.sd4x,col="red")
# points(rep(6,100),bootstrapped.sd6x,col="red")
# points(rep(8,100),bootstrapped.sd8x,col="red")
# points(rep(10,100),bootstrapped.sd10x,col="red")
# points(rep(12,100),bootstrapped.sd12x,col="red")



# plot(c(2,4,6,8,10,12),genome.means/genome.sd,xlim=c(1,14),pch=19,ylim=c(0,7))
# points(c(2,4,6,8,10,12),genome.means/genome.sd,type='l')
# points(rep(4,100),bootstrapped.mean4x/bootstrapped.sd4x,col="red")
# points(rep(6,100),bootstrapped.mean6x/bootstrapped.sd6x,col="red")
# points(rep(8,100),bootstrapped.mean8x/bootstrapped.sd8x,col="red")
# points(rep(10,100),bootstrapped.mean10x/bootstrapped.sd10x,col="red")
# points(rep(12,100),bootstrapped.mean12x/bootstrapped.sd12x,col="red")


qdmean4x<-sum((bootstrapped.mean4x-genome.means[2])^2)/B
qdmean6x<-sum((bootstrapped.mean6x-genome.means[3])^2)/B
qdmean8x<-sum((bootstrapped.mean8x-genome.means[4])^2)/B
qdmean10x<-sum((bootstrapped.mean10x-genome.means[5])^2)/B
qdmean12x<-sum((bootstrapped.mean12x-genome.means[6])^2)/B


qdmeanh1<-c(qdmean4x,qdmean6x,qdmean8x,qdmean10x,qdmean12x)
return(qdmeanh1)
}

qdmean<-retention.rate3(rr=1,B)
best.ateachploidy[3,]<-qdmean
	


points(c(4,6,8,10,12), best.ateachploidy[3,], cex=3,type="l",col=col1[3],lwd=2)





##############################################################################################################
##############Hypothesis 4: Thresholds monoploid number and Genome Size
retention.rate4<-function(rr,B){

bootstrapped.mean3x<-rep(0,B)
bootstrapped.sd3x<-rep(0,B)
bootstrapped.mean4x<-rep(0,B)
bootstrapped.sd4x<-rep(0,B)
bootstrapped.mean5x<-rep(0,B)
bootstrapped.sd5x<-rep(0,B)
bootstrapped.mean6x<-rep(0,B)
bootstrapped.sd6x<-rep(0,B)
bootstrapped.mean7x<-rep(0,B)
bootstrapped.sd7x<-rep(0,B)
bootstrapped.mean8x<-rep(0,B)
bootstrapped.sd8x<-rep(0,B)
bootstrapped.mean9x<-rep(0,B)
bootstrapped.sd9x<-rep(0,B)
bootstrapped.mean10x<-rep(0,B)
bootstrapped.sd10x<-rep(0,B)
bootstrapped.mean11x<-rep(0,B)
bootstrapped.sd11x<-rep(0,B)
bootstrapped.mean12x<-rep(0,B)
bootstrapped.sd12x<-rep(0,B)



#retention rate
for (i in 1:B){
	
	genome.sample<-genome.size[diploids]
	mono.sample<-mono.number[diploids]
	individualstoduplicate<-which((genome.sample< 53.4)&(mono.sample<31))	
		
	#Tetraploids
bootstrap.sample4x<-sample(x=genome.sample[individualstoduplicate], size=how.many4x, replace = TRUE, prob = NULL)*2*rr
bootstrapped.mean4x[i]<-mean(bootstrap.sample4x)
bootstrapped.sd4x[i]<-sd(bootstrap.sample4x)
	
	
	genome.sample<-genome.size[tetraploids]
	mono.sample<-mono.number[tetraploids]
	individualstoduplicate<-which((genome.sample<80.61)&(mono.sample<31))	
	
	#Hexaploids
bootstrap.sample6x<-sample(x=genome.sample[individualstoduplicate], size=how.many6x, replace = TRUE, prob = NULL)*1.5*rr
bootstrapped.mean6x[i]<-mean(bootstrap.sample6x)
bootstrapped.sd6x[i]<-sd(bootstrap.sample6x)


	genome.sample<-genome.size[hexaploids]
	mono.sample<-mono.number[hexaploids]
	individualstoduplicate<-which((genome.sample<78.71)&(mono.sample<31))
#Octoploids
bootstrap.sample8x<-sample(x=genome.sample[individualstoduplicate],size=how.many8x,replace=TRUE,prob=NULL)*1.33*rr
bootstrapped.mean8x[i]<-mean(bootstrap.sample8x)
bootstrapped.sd8x[i]<-sd(bootstrap.sample8x)


	genome.sample<-genome.size[octoploids]
	mono.sample<-mono.number[octoploids]
	individualstoduplicate<-which((genome.sample< 47.49)&(mono.sample<25))
	
	#Decaploids
bootstrap.sample10x<-sample(x=genome.sample[individualstoduplicate], size=how.many10x, replace = TRUE, prob = NULL)*1.25*rr
bootstrapped.mean10x[i]<-mean(bootstrap.sample10x)
bootstrapped.sd10x[i]<-sd(bootstrap.sample10x)

 
	genome.sample<-genome.size[decaploids]
	mono.sample<-mono.number[decaploids]
	individualstoduplicate<-which((genome.sample< 12.30)&(mono.sample<10))
	

bootstrap.sample12x<-sample(x=genome.sample[individualstoduplicate], size=how.many12x, replace = TRUE, prob = NULL)*1.2*rr
bootstrapped.mean12x[i]<-mean(bootstrap.sample12x)
bootstrapped.sd12x[i]<-sd(bootstrap.sample12x)
}

genome.means<-c(mean(genome.size[diploids]),mean(genome.size[tetraploids]),mean(genome.size[hexaploids]) ,mean(genome.size[octoploids]),mean(genome.size[decaploids]),mean(genome.size[duodecaploids]))

# plot(c(2,4,6,8,10,12),genome.means,xlim=c(1,14),pch=19,ylim=c(0,30))
# points(c(2,4,6,8,10,12),genome.means,type='l')
# points(rep(4,100),bootstrapped.mean4x,col="red")
# points(rep(6,100),bootstrapped.mean6x,col="red")
# points(rep(8,100),bootstrapped.mean8x,col="red")
# points(rep(10,100),bootstrapped.mean10x,col="red")
# points(rep(12,100),bootstrapped.mean12x,col="red")


genome.sd<-c(sd(genome.size[diploids]),sd(genome.size[tetraploids]),sd(genome.size[hexaploids]),sd(genome.size[octoploids]),sd(genome.size[decaploids]),sd(genome.size[duodecaploids]))

# plot(c(2,4,6,8,10,12),genome.sd,xlim=c(1,14),pch=19,ylim=c(0,30))
# points(c(2,4,6,8,10,12),genome.sd,type='l')
# points(rep(4,100),bootstrapped.sd4x,col="red")
# points(rep(6,100),bootstrapped.sd6x,col="red")
# points(rep(8,100),bootstrapped.sd8x,col="red")
# points(rep(10,100),bootstrapped.sd10x,col="red")
# points(rep(12,100),bootstrapped.sd12x,col="red")



# plot(c(2,4,6,8,10,12),genome.means/genome.sd,xlim=c(1,14),pch=19,ylim=c(0,7))
# points(c(2,4,6,8,10,12),genome.means/genome.sd,type='l')
# points(rep(4,100),bootstrapped.mean4x/bootstrapped.sd4x,col="red")
# points(rep(6,100),bootstrapped.mean6x/bootstrapped.sd6x,col="red")
# points(rep(8,100),bootstrapped.mean8x/bootstrapped.sd8x,col="red")
# points(rep(10,100),bootstrapped.mean10x/bootstrapped.sd10x,col="red")
# points(rep(12,100),bootstrapped.mean12x/bootstrapped.sd12x,col="red")


qdmean4x<-sum((bootstrapped.mean4x-genome.means[2])^2)/B
qdmean6x<-sum((bootstrapped.mean6x-genome.means[3])^2)/B
qdmean8x<-sum((bootstrapped.mean8x-genome.means[4])^2)/B
qdmean10x<-sum((bootstrapped.mean10x-genome.means[5])^2)/B
qdmean12x<-sum((bootstrapped.mean12x-genome.means[6])^2)/B


qdmeanh1<-c(qdmean4x,qdmean6x,qdmean8x,qdmean10x,qdmean12x)
return(qdmeanh1)
}

qdmean<-retention.rate4(rr=1,B)
best.ateachploidy[4,]<-qdmean
	


points(c(4,6,8,10,12), best.ateachploidy[4,], cex=3,type="l",col=col1[4],lwd=2)




##############################################################################################################
##############Hypothesis 5: Genome Size Threshold and Genome Downsizing
retention.rate5<-function(rr,B){

bootstrapped.mean3x<-rep(0,B)
bootstrapped.sd3x<-rep(0,B)
bootstrapped.mean4x<-rep(0,B)
bootstrapped.sd4x<-rep(0,B)
bootstrapped.mean5x<-rep(0,B)
bootstrapped.sd5x<-rep(0,B)
bootstrapped.mean6x<-rep(0,B)
bootstrapped.sd6x<-rep(0,B)
bootstrapped.mean7x<-rep(0,B)
bootstrapped.sd7x<-rep(0,B)
bootstrapped.mean8x<-rep(0,B)
bootstrapped.sd8x<-rep(0,B)
bootstrapped.mean9x<-rep(0,B)
bootstrapped.sd9x<-rep(0,B)
bootstrapped.mean10x<-rep(0,B)
bootstrapped.sd10x<-rep(0,B)
bootstrapped.mean11x<-rep(0,B)
bootstrapped.sd11x<-rep(0,B)
bootstrapped.mean12x<-rep(0,B)
bootstrapped.sd12x<-rep(0,B)



#retention rate
for (i in 1:B){
	
	genome.sample<-genome.size[diploids]
	mono.sample<-mono.number[diploids]
	individualstoduplicate<-which(genome.sample< 53.4)	
		
	#Tetraploids
bootstrap.sample4x<-sample(x=genome.sample[individualstoduplicate], size=how.many4x, replace = TRUE, prob = NULL)*2*rr
bootstrapped.mean4x[i]<-mean(bootstrap.sample4x)
bootstrapped.sd4x[i]<-sd(bootstrap.sample4x)
	
	
	genome.sample<-genome.size[tetraploids]
	mono.sample<-mono.number[tetraploids]
	individualstoduplicate<-which(genome.sample<80.61)
	
	#Hexaploids
bootstrap.sample6x<-sample(x=genome.sample[individualstoduplicate], size=how.many6x, replace = TRUE, prob = NULL)*1.5*rr
bootstrapped.mean6x[i]<-mean(bootstrap.sample6x)
bootstrapped.sd6x[i]<-sd(bootstrap.sample6x)


	genome.sample<-genome.size[hexaploids]
	mono.sample<-mono.number[hexaploids]
	individualstoduplicate<-which(genome.sample<78.71)
#Octoploids
bootstrap.sample8x<-sample(x=genome.sample[individualstoduplicate],size=how.many8x,replace=TRUE,prob=NULL)*1.33*rr
bootstrapped.mean8x[i]<-mean(bootstrap.sample8x)
bootstrapped.sd8x[i]<-sd(bootstrap.sample8x)


	genome.sample<-genome.size[octoploids]
	mono.sample<-mono.number[octoploids]
	individualstoduplicate<-which(genome.sample< 47.49)
	
	#Decaploids
bootstrap.sample10x<-sample(x=genome.sample[individualstoduplicate], size=how.many10x, replace = TRUE, prob = NULL)*1.25*rr
bootstrapped.mean10x[i]<-mean(bootstrap.sample10x)
bootstrapped.sd10x[i]<-sd(bootstrap.sample10x)

 
	genome.sample<-genome.size[decaploids]
	mono.sample<-mono.number[decaploids]
	individualstoduplicate<-which(genome.sample< 12.30)	

bootstrap.sample12x<-sample(x=genome.sample[individualstoduplicate], size=how.many12x, replace = TRUE, prob = NULL)*1.2*rr
bootstrapped.mean12x[i]<-mean(bootstrap.sample12x)
bootstrapped.sd12x[i]<-sd(bootstrap.sample12x)
}

genome.means<-c(mean(genome.size[diploids]),mean(genome.size[tetraploids]),mean(genome.size[hexaploids]) ,mean(genome.size[octoploids]),mean(genome.size[decaploids]),mean(genome.size[duodecaploids]))

# plot(c(2,4,6,8,10,12),genome.means,xlim=c(1,14),pch=19,ylim=c(0,30))
# points(c(2,4,6,8,10,12),genome.means,type='l')
# points(rep(4,100),bootstrapped.mean4x,col="red")
# points(rep(6,100),bootstrapped.mean6x,col="red")
# points(rep(8,100),bootstrapped.mean8x,col="red")
# points(rep(10,100),bootstrapped.mean10x,col="red")
# points(rep(12,100),bootstrapped.mean12x,col="red")


genome.sd<-c(sd(genome.size[diploids]),sd(genome.size[tetraploids]),sd(genome.size[hexaploids]),sd(genome.size[octoploids]),sd(genome.size[decaploids]),sd(genome.size[duodecaploids]))

# plot(c(2,4,6,8,10,12),genome.sd,xlim=c(1,14),pch=19,ylim=c(0,30))
# points(c(2,4,6,8,10,12),genome.sd,type='l')
# points(rep(4,100),bootstrapped.sd4x,col="red")
# points(rep(6,100),bootstrapped.sd6x,col="red")
# points(rep(8,100),bootstrapped.sd8x,col="red")
# points(rep(10,100),bootstrapped.sd10x,col="red")
# points(rep(12,100),bootstrapped.sd12x,col="red")



# plot(c(2,4,6,8,10,12),genome.means/genome.sd,xlim=c(1,14),pch=19,ylim=c(0,7))
# points(c(2,4,6,8,10,12),genome.means/genome.sd,type='l')
# points(rep(4,100),bootstrapped.mean4x/bootstrapped.sd4x,col="red")
# points(rep(6,100),bootstrapped.mean6x/bootstrapped.sd6x,col="red")
# points(rep(8,100),bootstrapped.mean8x/bootstrapped.sd8x,col="red")
# points(rep(10,100),bootstrapped.mean10x/bootstrapped.sd10x,col="red")
# points(rep(12,100),bootstrapped.mean12x/bootstrapped.sd12x,col="red")


qdmean4x<-sum((bootstrapped.mean4x-genome.means[2])^2)/B
qdmean6x<-sum((bootstrapped.mean6x-genome.means[3])^2)/B
qdmean8x<-sum((bootstrapped.mean8x-genome.means[4])^2)/B
qdmean10x<-sum((bootstrapped.mean10x-genome.means[5])^2)/B
qdmean12x<-sum((bootstrapped.mean12x-genome.means[6])^2)/B


qdmeanh1<-c(qdmean4x,qdmean6x,qdmean8x,qdmean10x,qdmean12x)
return(qdmeanh1)
}

min.qdmean5<-rep(0,long1) #Result of each hypothesis
bootstrapped.sim5<-rep(0,5)# Samples with each rate
for(j in 1:long1){
	qdmean<-retention.rate5(rr=ret.rates[j],B)
	bootstrapped.sim5<-rbind(bootstrapped.sim5,qdmean)
	min.qdmean5[j]<-sum(qdmean)	
}

write.table(bootstrapped.sim5,file="bootstrap5.txt",sep="\n",row.names=FALSE,col.names=FALSE)
#Plot the ploidy vs. MSE estimate 

#plot(c(2,4,6,8,10,12),rep(0,6),type="l",ylim=c(0,100), col="white",xlab="Ploidy",ylab="Mean Squared Error Estimate")
#for(j in 2:long1){
#	points(c(4,6,8,10,12), bootstrapped.sim1[j,], cex=3,type="l",col=col2[j])
#}
besthyp5<-which.min(min.qdmean5)
best.ratehyp5<-ret.rates[besthyp5] # 1- best.ratehyp 0.168
best.ateachploidy[5,]<-bootstrapped.sim5[besthyp5+1,]
points(c(4,6,8,10,12), best.ateachploidy[5,], cex=3,type="l",col=col1[5],lwd=2)



############
##############################################################################################################
##############Hypothesis 6: Thresholds monoploid number and Genome Size and Genome downsizing
retention.rate6<-function(rr,B){

bootstrapped.mean3x<-rep(0,B)
bootstrapped.sd3x<-rep(0,B)
bootstrapped.mean4x<-rep(0,B)
bootstrapped.sd4x<-rep(0,B)
bootstrapped.mean5x<-rep(0,B)
bootstrapped.sd5x<-rep(0,B)
bootstrapped.mean6x<-rep(0,B)
bootstrapped.sd6x<-rep(0,B)
bootstrapped.mean7x<-rep(0,B)
bootstrapped.sd7x<-rep(0,B)
bootstrapped.mean8x<-rep(0,B)
bootstrapped.sd8x<-rep(0,B)
bootstrapped.mean9x<-rep(0,B)
bootstrapped.sd9x<-rep(0,B)
bootstrapped.mean10x<-rep(0,B)
bootstrapped.sd10x<-rep(0,B)
bootstrapped.mean11x<-rep(0,B)
bootstrapped.sd11x<-rep(0,B)
bootstrapped.mean12x<-rep(0,B)
bootstrapped.sd12x<-rep(0,B)



#retention rate
for (i in 1:B){
	
	genome.sample<-genome.size[diploids]
	mono.sample<-mono.number[diploids]
	individualstoduplicate<-which((genome.sample< 53.4)&(mono.sample<31))	
		
	#Tetraploids
bootstrap.sample4x<-sample(x=genome.sample[individualstoduplicate], size=how.many4x, replace = TRUE, prob = NULL)*2*rr
bootstrapped.mean4x[i]<-mean(bootstrap.sample4x)
bootstrapped.sd4x[i]<-sd(bootstrap.sample4x)
	
	
	genome.sample<-genome.size[tetraploids]
	mono.sample<-mono.number[tetraploids]
	individualstoduplicate<-which((genome.sample<80.61)&(mono.sample<31))	
	
	#Hexaploids
bootstrap.sample6x<-sample(x=genome.sample[individualstoduplicate], size=how.many6x, replace = TRUE, prob = NULL)*1.5*rr
bootstrapped.mean6x[i]<-mean(bootstrap.sample6x)
bootstrapped.sd6x[i]<-sd(bootstrap.sample6x)


	genome.sample<-genome.size[hexaploids]
	mono.sample<-mono.number[hexaploids]
	individualstoduplicate<-which((genome.sample<78.71)&(mono.sample<31))
#Octoploids
bootstrap.sample8x<-sample(x=genome.sample[individualstoduplicate],size=how.many8x,replace=TRUE,prob=NULL)*1.33*rr
bootstrapped.mean8x[i]<-mean(bootstrap.sample8x)
bootstrapped.sd8x[i]<-sd(bootstrap.sample8x)


	genome.sample<-genome.size[octoploids]
	mono.sample<-mono.number[octoploids]
	individualstoduplicate<-which((genome.sample< 47.49)&(mono.sample<25))
	
	#Decaploids
bootstrap.sample10x<-sample(x=genome.sample[individualstoduplicate], size=how.many10x, replace = TRUE, prob = NULL)*1.25*rr
bootstrapped.mean10x[i]<-mean(bootstrap.sample10x)
bootstrapped.sd10x[i]<-sd(bootstrap.sample10x)

 
	genome.sample<-genome.size[decaploids]
	mono.sample<-mono.number[decaploids]
	individualstoduplicate<-which((genome.sample< 12.30)&(mono.sample<10))
	

bootstrap.sample12x<-sample(x=genome.sample[individualstoduplicate], size=how.many12x, replace = TRUE, prob = NULL)*1.2*rr
bootstrapped.mean12x[i]<-mean(bootstrap.sample12x)
bootstrapped.sd12x[i]<-sd(bootstrap.sample12x)
}

genome.means<-c(mean(genome.size[diploids]),mean(genome.size[tetraploids]),mean(genome.size[hexaploids]) ,mean(genome.size[octoploids]),mean(genome.size[decaploids]),mean(genome.size[duodecaploids]))

# plot(c(2,4,6,8,10,12),genome.means,xlim=c(1,14),pch=19,ylim=c(0,30))
# points(c(2,4,6,8,10,12),genome.means,type='l')
# points(rep(4,100),bootstrapped.mean4x,col="red")
# points(rep(6,100),bootstrapped.mean6x,col="red")
# points(rep(8,100),bootstrapped.mean8x,col="red")
# points(rep(10,100),bootstrapped.mean10x,col="red")
# points(rep(12,100),bootstrapped.mean12x,col="red")


genome.sd<-c(sd(genome.size[diploids]),sd(genome.size[tetraploids]),sd(genome.size[hexaploids]),sd(genome.size[octoploids]),sd(genome.size[decaploids]),sd(genome.size[duodecaploids]))

# plot(c(2,4,6,8,10,12),genome.sd,xlim=c(1,14),pch=19,ylim=c(0,30))
# points(c(2,4,6,8,10,12),genome.sd,type='l')
# points(rep(4,100),bootstrapped.sd4x,col="red")
# points(rep(6,100),bootstrapped.sd6x,col="red")
# points(rep(8,100),bootstrapped.sd8x,col="red")
# points(rep(10,100),bootstrapped.sd10x,col="red")
# points(rep(12,100),bootstrapped.sd12x,col="red")



# plot(c(2,4,6,8,10,12),genome.means/genome.sd,xlim=c(1,14),pch=19,ylim=c(0,7))
# points(c(2,4,6,8,10,12),genome.means/genome.sd,type='l')
# points(rep(4,100),bootstrapped.mean4x/bootstrapped.sd4x,col="red")
# points(rep(6,100),bootstrapped.mean6x/bootstrapped.sd6x,col="red")
# points(rep(8,100),bootstrapped.mean8x/bootstrapped.sd8x,col="red")
# points(rep(10,100),bootstrapped.mean10x/bootstrapped.sd10x,col="red")
# points(rep(12,100),bootstrapped.mean12x/bootstrapped.sd12x,col="red")


qdmean4x<-sum((bootstrapped.mean4x-genome.means[2])^2)/B
qdmean6x<-sum((bootstrapped.mean6x-genome.means[3])^2)/B
qdmean8x<-sum((bootstrapped.mean8x-genome.means[4])^2)/B
qdmean10x<-sum((bootstrapped.mean10x-genome.means[5])^2)/B
qdmean12x<-sum((bootstrapped.mean12x-genome.means[6])^2)/B


qdmeanh1<-c(qdmean4x,qdmean6x,qdmean8x,qdmean10x,qdmean12x)
return(qdmeanh1)
}

min.qdmean6<-rep(0,long1) #Result of each hypothesis
bootstrapped.sim6<-rep(0,5)# Samples with each rate
for(j in 1:long1){
	qdmean<-retention.rate6(rr=ret.rates[j],B)
	bootstrapped.sim6<-rbind(bootstrapped.sim6,qdmean)
	min.qdmean6[j]<-sum(qdmean)	
}

write.table(bootstrapped.sim6,file="bootstrap6.txt",sep="\n",row.names=FALSE,col.names=FALSE)
#Plot the ploidy vs. MSE estimate 

#plot(c(2,4,6,8,10,12),rep(0,6),type="l",ylim=c(0,100), col="white",xlab="Ploidy",ylab="Mean Squared Error Estimate")
#for(j in 2:long1){
#	points(c(4,6,8,10,12), bootstrapped.sim1[j,], cex=3,type="l",col=col2[j])
#}
besthyp6<-which.min(min.qdmean6)
best.ratehyp6<-ret.rates[besthyp6] # 1- best.ratehyp 0.168
best.ateachploidy[6,]<-bootstrapped.sim6[besthyp6+1,]
points(c(4,6,8,10,12), best.ateachploidy[6,], cex=3,type="l",col=col1[6],lwd=2)
legend(2,38,legend=c("GD","MN Threshold","GS Threshold","MN+GS Threshold","GS Threshold+GD","MN+GS Threshold+GD"),pt.cex=0.5,col=col1,lwd=3)

 write.table(best.ateachploidy,file="bestateachploidy.txt",sep="\n",row.names=FALSE,col.names=FALSE)
 
########################################################################### 
 ####################PLOTS
 
 plot(c(2,4,6,8,10,12),rep(0,6),type="l",ylim=c(0,40), col="white",xlab="Ploidy",ylab="Mean Squared Error Estimate")
 for(i in 1:6){
 points(c(4,6,8,10,12), best.ateachploidy[i,], cex=3,type="l",col=col1[i],lwd=2)
}
legend(2,38,legend=c("GD","MN Threshold","GS Threshold","MN+GS Threshold","GS Threshold+GD","MN+GS Threshold+GD"),pt.cex=0.5,col=col1,lwd=3)


  plot(c(2,4,6,8,10,12),rep(0,6),type="l",ylim=c(0,80), col="white",xlab="Ploidy",ylab="Mean Squared Error Estimate")
for(j in 2:long1+1){
	points(c(4,6,8,10,12), bootstrapped.sim1[j,], cex=3,type="l",col=col2[j])
}
points(c(4,6,8,10,12), best.ateachploidy[1,], cex=3,type="l",col=col1[1],lwd=3)
}
legend(2,75,legend=c("Best GD rate 0.301"),pt.cex=0.5,col=col1[1],lwd=3)


#######################
plot(c(2,4,6,8,10,12),rep(0,6),type="l",ylim=c(0,80), col="white",xlab="Ploidy",ylab="Mean Squared Error Estimate")

for(j in 2:long1+1){
	points(c(4,6,8,10,12), bootstrapped.sim5[j,], cex=3,type="l",col=col2[j])
}
points(c(4,6,8,10,12), best.ateachploidy[5,], cex=3,type="l",col=col1[5],lwd=3)
}
legend(2,75,legend=c("Best GD rate + GS Threshold 0.168"),pt.cex=0.5,col=col1[5],lwd=3)

#######################
plot(c(2,4,6,8,10,12),rep(0,6),type="l",ylim=c(0,80), col="white",xlab="Ploidy",ylab="Mean Squared Error Estimate")

for(j in 2:long1+1){
	points(c(4,6,8,10,12), bootstrapped.sim6[j,], cex=3,type="l",col=col2[j])
}
points(c(4,6,8,10,12), best.ateachploidy[6,], cex=3,type="l",col=col1[6],lwd=3)
}
legend(2,75,legend=c("Best GD rate + GS and MN Threshold 0.154"),pt.cex=0.5,col=col1[6],lwd=3)

#####
plot(1-ret.rates,apply(bootstrapped.sim1[-1,],1,sum),col=col1[1],type="l",xlab="Genome Downsizing Rates",ylab="Total MSE",lwd=1.2)

points(1-best.ratehyp1,sum(best.ateachploidy[1,]),col=col1[1],pch=10)

points(1-ret.rates,apply(bootstrapped.sim5[-1,],1,sum),col=col1[5],type="l",lwd=1.2)
points(1-best.ratehyp5,sum(best.ateachploidy[5,]),col=col1[5],pch=10)
points(1-ret.rates,apply(bootstrapped.sim6[-1,],1,sum),col=col1[6],type="l",lwd=1.2)
points(1-best.ratehyp6,sum(best.ateachploidy[6,]),col=col1[6],pch=10)

legend(0.1,255, legend=c("GD 0.301","GD + GS Threshold 0.168","GD + GS and MN Threshold 0.154"),pt.cex=0.5, col=col1[c(1,5,6)],lwd=3)

