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
non.repeated<-which(angiocvalues$Minimum==1)#7124
polyploidseries<-angiocvalues[,15]

# 1-Family, 2-Genus, 3-Species, 4-Authority, 5-Chromosome.number, 6-Ploidy.level, 7-Estimation.method, 8-Voucher, 9-X1C..pg., 10-Original.Reference,
# 11-Lifecycle, 12-Group, 13-Paper, 14-Name, 15-History.of.Polyploidy, 16-Repeated.Sample, 17-Non.multiple, 18-Minimum, 19-Missing.Data, 20-Possible.Hybrid, 
#21-Both, 22-Chromploid
 

genome.size<-angiocvalues[non.repeated,9]
ploidy.genome<-angiocvalues[non.repeated,6]
mono.number<-angiocvalues[non.repeated,5]/angiocvalues[non.repeated,6]
genome.boxplot<-boxplot(genome.size~ploidy.genome)
genus.genome<-as.factor(angiocvalues[non.repeated,2]) #1631
family.genome<-as.factor(angiocvalues[non.repeated,1])#248
how.manygenus<-length(levels(genus.genome))
minitable<-angiocvalues[non.repeated,c(2,6,9)]
genusbyploidy=table(minitable[,1],minitable[,2])
#write.table(genusbyploidy,file="genusbyploidy.txt",sep=",",row.names=TRUE, col.names=TRUE)

genus.names<-levels(genus.genome)
dim(genusbyploidy)
apply(genusbyploidy,2,FUN="sum")

diploid.list<-list()
for (i in 1:how.manygenus){
	aux1<- which((genus.genome==genus.names[i])&(ploidy.genome==2))
	specific.genome<-genome.size[aux1]
	if(length(aux1)>0){
		diploid.list<-c(diploid.list,list(specific.genome))
	}
}
how.many2xbygenus<-length(diploid.list)
averagebygenus.2x<-rep(0,how.many2xbygenus)
for(j in 1:how.many2xbygenus){
	averagebygenus.2x[j]<-mean(diploid.list[[j]])
}
average2x<-mean(averagebygenus.2x)

triploid.list<-list()
for (i in 1:how.manygenus){
	aux1<- which((genus.genome==genus.names[i])&(ploidy.genome==3))
	specific.genome<-genome.size[aux1]
	if(length(aux1)>0){
		triploid.list<-c(triploid.list,list(specific.genome))
	}
}
how.many3xbygenus<-length(triploid.list)
averagebygenus.3x<-rep(0,how.many3xbygenus)
for(j in 1:how.many3xbygenus){
	averagebygenus.3x[j]<-mean(triploid.list[[j]])
}
average3x<-mean(averagebygenus.3x)



tetraploid.list<-list()
for (i in 1:how.manygenus){
	aux1<- which((genus.genome==genus.names[i])&(ploidy.genome==4))
	specific.genome<-genome.size[aux1]
	if(length(aux1)>0){
		tetraploid.list<-c(tetraploid.list,list(specific.genome))
	}
}
how.many4xbygenus<-length(tetraploid.list)
averagebygenus.4x<-rep(0,how.many4xbygenus)
for(j in 1:how.many4xbygenus){
	averagebygenus.4x[j]<-mean(tetraploid.list[[j]])
}
average4x<-mean(averagebygenus.4x)





pentaploid.list<-list()
for (i in 1:how.manygenus){
	aux1<- which((genus.genome==genus.names[i])&(ploidy.genome==5))
	specific.genome<-genome.size[aux1]
	if(length(aux1)>0){
		pentaploid.list<-c(pentaploid.list,list(specific.genome))
	}
}
how.many5xbygenus<-length(pentaploid.list)

averagebygenus.5x<-rep(0,how.many5xbygenus)
for(j in 1:how.many5xbygenus){
	averagebygenus.5x[j]<-mean(pentaploid.list[[j]])
}
average5x<-mean(averagebygenus.5x)




hexaploid.list<-list()
for (i in 1:how.manygenus){
	aux1<- which((genus.genome==genus.names[i])&(ploidy.genome==6))
	specific.genome<-genome.size[aux1]
	if(length(aux1)>0){
		hexaploid.list<-c(hexaploid.list,list(specific.genome))
	}
}
how.many6xbygenus<-length(hexaploid.list)
averagebygenus.6x<-rep(0,how.many6xbygenus)
for(j in 1:how.many6xbygenus){
	averagebygenus.6x[j]<-mean(hexaploid.list[[j]])
}
average6x<-mean(averagebygenus.6x)



heptaploid.list<-list()
for (i in 1:how.manygenus){
	aux1<- which((genus.genome==genus.names[i])&(ploidy.genome==7))
	specific.genome<-genome.size[aux1]
	if(length(aux1)>0){
		heptaploid.list<-c(heptaploid.list,list(specific.genome))
	}
}
how.many7xbygenus<-length(heptaploid.list)
averagebygenus.7x<-rep(0,how.many7xbygenus)
for(j in 1:how.many7xbygenus){
	averagebygenus.7x[j]<-mean(heptaploid.list[[j]])
}
average7x<-mean(averagebygenus.7x)


octoploid.list<-list()
for (i in 1:how.manygenus){
	aux1<- which((genus.genome==genus.names[i])&(ploidy.genome==8))
	specific.genome<-genome.size[aux1]
	if(length(aux1)>0){
		octoploid.list<-c(octoploid.list,list(specific.genome))
	}
}
how.many8xbygenus<-length(octoploid.list)
averagebygenus.8x<-rep(0,how.many8xbygenus)
for(j in 1:how.many8xbygenus){
	averagebygenus.8x[j]<-mean(octoploid.list[[j]])
}
average8x<-mean(averagebygenus.8x)




enneaploid.list<-list()
for (i in 1:how.manygenus){
	aux1<- which((genus.genome==genus.names[i])&(ploidy.genome==9))
	specific.genome<-genome.size[aux1]
	if(length(aux1)>0){
		enneaploid.list<-c(enneaploid.list,list(specific.genome))
	}
}
how.many9xbygenus<-length(enneaploid.list)
averagebygenus.9x<-rep(0,how.many9xbygenus)
for(j in 1:how.many9xbygenus){
	averagebygenus.9x[j]<-mean(enneaploid.list[[j]])
}
average9x<-mean(averagebygenus.9x)




decaploid.list<-list()
for (i in 1:how.manygenus){
	aux1<- which((genus.genome==genus.names[i])&(ploidy.genome==10))
	specific.genome<-genome.size[aux1]
	if(length(aux1)>0){
		decaploid.list<-c(decaploid.list,list(specific.genome))
	}
}
how.many10xbygenus<-length(decaploid.list)
averagebygenus.10x<-rep(0,how.many10xbygenus)
for(j in 1:how.many10xbygenus){
	averagebygenus.10x[j]<-mean(decaploid.list[[j]])
}
average10x<-mean(averagebygenus.10x)


undecaploid.list<-list()
for (i in 1:how.manygenus){
	aux1<- which((genus.genome==genus.names[i])&(ploidy.genome==11))
	specific.genome<-genome.size[aux1]
	if(length(aux1)>0){
		undecaploid.list<-c(undecaploid.list,list(specific.genome))
	}
}
how.many11xbygenus<-length(undecaploid.list)
averagebygenus.11x<-rep(0,how.many11xbygenus)
for(j in 1:how.many11xbygenus){
	averagebygenus.11x[j]<-mean(undecaploid.list[[j]])
}
average11x<-mean(averagebygenus.11x)

duodecaploid.list<-list()
for (i in 1:how.manygenus){
	aux1<- which((genus.genome==genus.names[i])&(ploidy.genome==12))
	specific.genome<-genome.size[aux1]
	if(length(aux1)>0){
		duodecaploid.list<-c(duodecaploid.list,list(specific.genome))
	}
}
how.many12xbygenus<-length(duodecaploid.list)
averagebygenus.12x<-rep(0,how.many12xbygenus)
for(j in 1:how.many12xbygenus){
	averagebygenus.12x[j]<-mean(duodecaploid.list[[j]])
}
average12x<-mean(averagebygenus.12x)

averagesbygenus<-c(average2x,average3x,average4x,average5x,average6x,average7x,average8x,average9x,average10x,average11x,average12x)
points(1:11, averagesbygenus, type="l", lwd=2, col="orange")
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
# ecdf(genome.size[diploids])(53.4)
# ecdf(genome.size[triploids])(77.55)
# ecdf(genome.size[tetraploids])(80.6)
# ecdf(genome.size[pentaploids])(26.65)
# ecdf(genome.size[hexaploids])(78.7)
# ecdf(genome.size[heptaploids])(6.93)
# ecdf(genome.size[octoploids])(47.48)
# ecdf(genome.size[enneaploids])(20.28)
# ecdf(genome.size[decaploids])(12.29)
# ecdf(genome.size[undecaploids])(1.12)
# ecdf(genome.size[duodecaploids])(12.84)

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

#### Sample sizes by genus
sample.size2x<-min(how.many2xbygenus,how.many4x) #893 but we needed 914 for 4x 
sample.size4x<-min(how.many4xbygenus,how.many6x) #287 correct number of 6x
sample.size6x<-min(how.many6xbygenus,how.many8x) #99 correct number of 8x
sample.size8x<-min(how.many8xbygenus,how.many10x) #23 correct number of 10x
sample.size10x<-min(how.many10xbygenus,how.many12x)#13 but we needed 19 for 12x

 #retention rate
for (i in 1:B){
	
	#Creating new 4x from sample 2x by genus
	
	      diploid.sample<-rep(0,sample.size2x) #for 4x all genera of 2x need to be sampled
	for(j in 1:sample.size2x){
		aux1<-sample(diploid.list[[j]],1,replace=FALSE)
		diploid.sample[j]<-aux1
	}
		
	
bootstrap.sample4x<-diploid.sample*2*rr
bootstrapped.mean4x[i]<-mean(bootstrap.sample4x) # Mean of the new bootstrap sample
bootstrapped.sd4x[i]<-sd(bootstrap.sample4x) # Standard error of the new bootstrap sample
	
	
	# Creating new 6x from sample 4x by genus
	
	sample.genus<-sort(sample(1:how.many4xbygenus,sample.size4x,replace=FALSE)) #for 6x only some 4x genera need to be sampled
	tetraploid.sample<-rep(0,sample.size4x)
	for(j in 1:sample.size4x){
		aux1<-sample(tetraploid.list[[sample.genus[j]]],1,replace=FALSE)
		tetraploid.sample[j]<-aux1
	}
	
bootstrap.sample6x<-tetraploid.sample*1.5*rr
bootstrapped.mean6x[i]<-mean(bootstrap.sample6x)
bootstrapped.sd6x[i]<-sd(bootstrap.sample6x)


# Creating new 8x from sample 6x by genus
	
	sample.genus<-sort(sample(1:how.many6xbygenus,sample.size6x,replace=FALSE)) #for 8x only some 6x genera need to be sampled
	hexaploid.sample<-rep(0,sample.size6x)
	for(j in 1:sample.size6x){
		aux1<-sample(hexaploid.list[[sample.genus[j]]],1,replace=FALSE)
		hexaploid.sample[j]<-aux1
	}
	

bootstrap.sample8x<-hexaploid.sample*1.33*rr
bootstrapped.mean8x[i]<-mean(bootstrap.sample8x)
bootstrapped.sd8x[i]<-sd(bootstrap.sample8x)


	# Creating new 10x from sample 8x by genus
	
	sample.genus<-sort(sample(1:how.many8xbygenus,sample.size8x,replace=FALSE)) #for 10x only some 8x genera need to be sampled
	octoploid.sample<-rep(0,sample.size8x)
	for(j in 1:sample.size8x){
		aux1<-sample(octoploid.list[[sample.genus[j]]],1,replace=FALSE)
		octoploid.sample[j]<-aux1
	}
	
bootstrap.sample10x<-octoploid.sample*1.25*rr
bootstrapped.mean10x[i]<-mean(bootstrap.sample10x)
bootstrapped.sd10x[i]<-sd(bootstrap.sample10x)

 
	# Creating new 12x from sample 10x by genus
	
	decaploid.sample<-rep(0,sample.size10x)# for 12x all of the genera of 10x need to be sampled because they are not enough
	for(j in 1:sample.size10x){
		aux1<-sample(decaploid.list[[j]],1,replace=FALSE)
		decaploid.sample[j]<-aux1
	}
	

	

bootstrap.sample12x<-decaploid.sample*1.2*rr
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

write.table(bootstrapped.sim1,file="bootstrapbygenusgd1.txt",sep="\n",row.names=FALSE,col.names=FALSE)
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

#### I need to create the sample first that doesn't pass that monoploid number
diploid.listmn<-list()
for (i in 1:how.manygenus){
	aux1<- which((genus.genome==genus.names[i])&(ploidy.genome==2)&(mono.number<31))
	specific.genome<-genome.size[aux1]
	if(length(aux1)>0){
		diploid.listmn<-c(diploid.listmn,list(specific.genome))
	}
}
how.many2xbygenusmn<-length(diploid.listmn)



tetraploid.listmn<-list()
for (i in 1:how.manygenus){
	aux1<- which((genus.genome==genus.names[i])&(ploidy.genome==4)&(mono.number<31))
	specific.genome<-genome.size[aux1]
	if(length(aux1)>0){
		tetraploid.listmn<-c(tetraploid.listmn,list(specific.genome))
	}
}
how.many4xbygenusmn<-length(tetraploid.listmn)


hexaploid.listmn<-list()
for (i in 1:how.manygenus){
	aux1<- which((genus.genome==genus.names[i])&(ploidy.genome==6)&(mono.number<31))
	specific.genome<-genome.size[aux1]
	if(length(aux1)>0){
		hexaploid.listmn<-c(hexaploid.listmn,list(specific.genome))
	}
}
how.many6xbygenusmn<-length(hexaploid.listmn)


octoploid.listmn<-list()
for (i in 1:how.manygenus){
	aux1<- which((genus.genome==genus.names[i])&(ploidy.genome==8)&(mono.number<25))
	specific.genome<-genome.size[aux1]
	if(length(aux1)>0){
		octoploid.listmn<-c(octoploid.listmn,list(specific.genome))
	}
}
how.many8xbygenusmn<-length(octoploid.listmn)

decaploid.listmn<-list()
for (i in 1:how.manygenus){
	aux1<- which((genus.genome==genus.names[i])&(ploidy.genome==10)&(mono.number<10))
	specific.genome<-genome.size[aux1]
	if(length(aux1)>0){
		decaploid.listmn<-c(decaploid.listmn,list(specific.genome))
	}
}
how.many10xbygenusmn<-length(decaploid.listmn)


#### Sample sizes by genus
sample.size2xmn<-min(how.many2xbygenusmn,how.many4x) #855 but we needed 914 for 4x 
sample.size4xmn<-min(how.many4xbygenusmn,how.many6x) #287 correct number of 6x
sample.size6xmn<-min(how.many6xbygenusmn,how.many8x) #99 correct number of 8x
sample.size8xmn<-min(how.many8xbygenusmn,how.many10x) #27 correct number of 10x
sample.size10xmn<-min(how.many10xbygenusmn,how.many12x)#12 but we needed 19 for 12x


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
	#### Creating new tetraploids
	diploid.sample<-rep(0,sample.size2xmn) #for 4x all genera of 2x need to be sampled
	for(j in 1:sample.size2xmn){
		aux1<-sample(diploid.listmn[[j]],1,replace=FALSE)
		diploid.sample[j]<-aux1
	}
	bootstrap.sample4x<-diploid.sample*2*rr
	bootstrapped.mean4x[i]<-mean(bootstrap.sample4x)
	bootstrapped.sd4x[i]<-sd(bootstrap.sample4x)
	
	###Creating new hexaploids
	sample.genus<-sort(sample(1:how.many4xbygenusmn,sample.size4xmn,replace=FALSE)) #for 6x only some 4x genera need to be sampled
	tetraploid.sample<-rep(0,sample.size4xmn)
	for(j in 1:sample.size4xmn){
		aux1<-sample(tetraploid.listmn[[sample.genus[j]]],1,replace=FALSE)
		tetraploid.sample[j]<-aux1
	}
	bootstrap.sample6x<-tetraploid.sample*1.5*rr
	bootstrapped.mean6x[i]<-mean(bootstrap.sample6x)
	bootstrapped.sd6x[i]<-sd(bootstrap.sample6x)
	
	###Creating new octoploids
	sample.genus<-sort(sample(1:how.many6xbygenusmn,sample.size6xmn,replace=FALSE)) #for 8x only some 6x genera need to be sampled
	hexaploid.sample<-rep(0,sample.size6xmn)
	for(j in 1:sample.size6xmn){
		aux1<-sample(hexaploid.listmn[[sample.genus[j]]],1,replace=FALSE)
		hexaploid.sample[j]<-aux1
	}
	
bootstrap.sample8x<-hexaploid.sample*1.33*rr
bootstrapped.mean8x[i]<-mean(bootstrap.sample8x)
bootstrapped.sd8x[i]<-sd(bootstrap.sample8x)
	
######## Creating new decaploids
	sample.genus<-sort(sample(1:how.many8xbygenusmn,sample.size8xmn,replace=FALSE)) #for 8x only some 6x genera need to be sampled
	octoploid.sample<-rep(0,sample.size8xmn)
	for(j in 1:sample.size8xmn){
		aux1<-sample(octoploid.listmn[[sample.genus[j]]],1,replace=FALSE)
		octoploid.sample[j]<-aux1
	}
	
bootstrap.sample10x<-octoploid.sample*1.25*rr
bootstrapped.mean10x[i]<-mean(bootstrap.sample10x)
bootstrapped.sd10x[i]<-sd(bootstrap.sample10x)
	

###### Creating new duodecaploids	
	decaploid.sample<-rep(0,sample.size10xmn) #for 4x all genera of 2x need to be sampled
	for(j in 1:sample.size10xmn){
		aux1<-sample(decaploid.listmn[[j]],1,replace=FALSE)
		decaploid.sample[j]<-aux1
	}

	 bootstrap.sample12x<-decaploid*1.2*rr
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
#### I need to create the sample first that doesn't pass that genome size
diploid.listgs<-list()
for (i in 1:how.manygenus){
	aux1<- which((genus.genome==genus.names[i])&(ploidy.genome==2)&(genome.size<53.4))
	specific.genome<-genome.size[aux1]
	if(length(aux1)>0){
		diploid.listgs<-c(diploid.listgs,list(specific.genome))
	}
}
how.many2xbygenusgs<-length(diploid.listgs)



tetraploid.listgs<-list()
for (i in 1:how.manygenus){
	aux1<- which((genus.genome==genus.names[i])&(ploidy.genome==4)&(genome.size<80.61))
	specific.genome<-genome.size[aux1]
	if(length(aux1)>0){
		tetraploid.listgs<-c(tetraploid.listgs,list(specific.genome))
	}
}
how.many4xbygenusgs<-length(tetraploid.listgs)


hexaploid.listgs<-list()
for (i in 1:how.manygenus){
	aux1<- which((genus.genome==genus.names[i])&(ploidy.genome==6)&(genome.size<78.71))
	specific.genome<-genome.size[aux1]
	if(length(aux1)>0){
		hexaploid.listgs<-c(hexaploid.listgs,list(specific.genome))
	}
}
how.many6xbygenusgs<-length(hexaploid.listgs)


octoploid.listgs<-list()
for (i in 1:how.manygenus){
	aux1<- which((genus.genome==genus.names[i])&(ploidy.genome==8)&(genome.size<47.49))
	specific.genome<-genome.size[aux1]
	if(length(aux1)>0){
		octoploid.listgs<-c(octoploid.listgs,list(specific.genome))
	}
}
how.many8xbygenusgs<-length(octoploid.listgs)

decaploid.listgs<-list()
for (i in 1:how.manygenus){
	aux1<- which((genus.genome==genus.names[i])&(ploidy.genome==10)&(genome.size<12.30))
	specific.genome<-genome.size[aux1]
	if(length(aux1)>0){
		decaploid.listgs<-c(decaploid.listgs,list(specific.genome))
	}
}
how.many10xbygenusgs<-length(decaploid.listgs)


#### Sample sizes by genus
sample.size2xgs<-min(how.many2xbygenusgs,how.many4x) #891 but we needed 914 for 4x 
sample.size4xgs<-min(how.many4xbygenusgs,how.many6x) #287 correct number of 6x
sample.size6xgs<-min(how.many6xbygenusgs,how.many8x) #99 correct number of 8x
sample.size8xgs<-min(how.many8xbygenusgs,how.many10x) #27 correct number of 10x
sample.size10xgs<-min(how.many10xbygenusgs,how.many12x)#10 but we needed 19 for 12x


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
	#### Creating new tetraploids
	diploid.sample<-rep(0,sample.size2xgs) #for 4x all genera of 2x need to be sampled
	for(j in 1:sample.size2xgs){
		aux1<-sample(diploid.listgs[[j]],1,replace=FALSE)
		diploid.sample[j]<-aux1
	}
	bootstrap.sample4x<-diploid.sample*2*rr
	bootstrapped.mean4x[i]<-mean(bootstrap.sample4x)
	bootstrapped.sd4x[i]<-sd(bootstrap.sample4x)
	
	###Creating new hexaploids
	sample.genus<-sort(sample(1:how.many4xbygenusgs,sample.size4xgs,replace=FALSE)) #for 6x only some 4x genera need to be sampled
	tetraploid.sample<-rep(0,sample.size4xgs)
	for(j in 1:sample.size4xgs){
		aux1<-sample(tetraploid.listgs[[sample.genus[j]]],1,replace=FALSE)
		tetraploid.sample[j]<-aux1
	}
	bootstrap.sample6x<-tetraploid.sample*1.5*rr
	bootstrapped.mean6x[i]<-mean(bootstrap.sample6x)
	bootstrapped.sd6x[i]<-sd(bootstrap.sample6x)
	
	###Creating new octoploids
	sample.genus<-sort(sample(1:how.many6xbygenusgs,sample.size6xgs,replace=FALSE)) #for 8x only some 6x genera need to be sampled
	hexaploid.sample<-rep(0,sample.size6xgs)
	for(j in 1:sample.size6xgs){
		aux1<-sample(hexaploid.listgs[[sample.genus[j]]],1,replace=FALSE)
		hexaploid.sample[j]<-aux1
	}
	
bootstrap.sample8x<-hexaploid.sample*1.33*rr
bootstrapped.mean8x[i]<-mean(bootstrap.sample8x)
bootstrapped.sd8x[i]<-sd(bootstrap.sample8x)
	
######## Creating new decaploids
	sample.genus<-sort(sample(1:how.many8xbygenusgs,sample.size8xgs,replace=FALSE)) #for 8x only some 6x genera need to be sampled
	octoploid.sample<-rep(0,sample.size8xgs)
	for(j in 1:sample.size8xgs){
		aux1<-sample(octoploid.listgs[[sample.genus[j]]],1,replace=FALSE)
		octoploid.sample[j]<-aux1
	}
	
bootstrap.sample10x<-octoploid.sample*1.25*rr
bootstrapped.mean10x[i]<-mean(bootstrap.sample10x)
bootstrapped.sd10x[i]<-sd(bootstrap.sample10x)
	

###### Creating new duodecaploids	
	decaploid.sample<-rep(0,sample.size10xgs) #for 4x all genera of 2x need to be sampled
	for(j in 1:sample.size10xgs){
		aux1<-sample(decaploid.listgs[[j]],1,replace=FALSE)
		decaploid.sample[j]<-aux1
	}

	 bootstrap.sample12x<-decaploid*1.2*rr
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

qdmean<-retention.rate3(rr=1,B)
best.ateachploidy[3,]<-qdmean
	


points(c(4,6,8,10,12), best.ateachploidy[3,], cex=3,type="l",col=col1[3],lwd=2)





##############################################################################################################
##############Hypothesis 4: Thresholds monoploid number and Genome Size
#### I need to create the sample first that doesn't pass that genome size and monoploid number
diploid.listgsmn<-list()
for (i in 1:how.manygenus){
	aux1<- which((genus.genome==genus.names[i])&(ploidy.genome==2)&(genome.size<53.4)&(mono.number<31))
	specific.genome<-genome.size[aux1]
	if(length(aux1)>0){
		diploid.listgsmn<-c(diploid.listgsmn,list(specific.genome))
	}
}
how.many2xbygenusgsmn<-length(diploid.listgsmn)

tetraploid.listgsmn<-list()
for (i in 1:how.manygenus){
	aux1<- which((genus.genome==genus.names[i])&(ploidy.genome==4)&(genome.size<80.61)&(mono.number<31))
	specific.genome<-genome.size[aux1]
	if(length(aux1)>0){
		tetraploid.listgsmn<-c(tetraploid.listgsmn,list(specific.genome))
	}
}
how.many4xbygenusgsmn<-length(tetraploid.listgsmn)


hexaploid.listgsmn<-list()
for (i in 1:how.manygenus){
	aux1<- which((genus.genome==genus.names[i])&(ploidy.genome==6)&(genome.size<78.71)&(mono.number<31))
	specific.genome<-genome.size[aux1]
	if(length(aux1)>0){
		hexaploid.listgsmn<-c(hexaploid.listgsmn,list(specific.genome))
	}
}
how.many6xbygenusgsmn<-length(hexaploid.listgsmn)


octoploid.listgsmn<-list()
for (i in 1:how.manygenus){
	aux1<- which((genus.genome==genus.names[i])&(ploidy.genome==8)&(genome.size<47.49)&(mono.number<25))
	specific.genome<-genome.size[aux1]
	if(length(aux1)>0){
		octoploid.listgsmn<-c(octoploid.listgsmn,list(specific.genome))
	}
}
how.many8xbygenusgsmn<-length(octoploid.listgsmn)

decaploid.listgsmn<-list()
for (i in 1:how.manygenus){
	aux1<- which((genus.genome==genus.names[i])&(ploidy.genome==10)&(genome.size<12.30))
	specific.genome<-genome.size[aux1]
	if(length(aux1)>0){
		decaploid.listgsmn<-c(decaploid.listgsmn,list(specific.genome))
	}
}
how.many10xbygenusgsmn<-length(decaploid.listgsmn)


#### Sample sizes by genus
sample.size2xgsmn<-min(how.many2xbygenusgsmn,how.many4x) #853 but we needed 914 for 4x 
sample.size4xgsmn<-min(how.many4xbygenusgsmn,how.many6x) #287 correct number of 6x
sample.size6xgsmn<-min(how.many6xbygenusgsmn,how.many8x) #99 correct number of 8x
sample.size8xgsmn<-min(how.many8xbygenusgsmn,how.many10x) #27 correct number of 10x
sample.size10xgsmn<-min(how.many10xbygenusgsmn,how.many12x)#10 but we needed 19 for 12x


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
	#### Creating new tetraploids
	diploid.sample<-rep(0,sample.size2xgsmn) #for 4x all genera of 2x need to be sampled
	for(j in 1:sample.size2xgsmn){
		aux1<-sample(diploid.listgs[[j]],1,replace=FALSE)
		diploid.sample[j]<-aux1
	}
	bootstrap.sample4x<-diploid.sample*2*rr
	bootstrapped.mean4x[i]<-mean(bootstrap.sample4x)
	bootstrapped.sd4x[i]<-sd(bootstrap.sample4x)
	
	###Creating new hexaploids
	sample.genus<-sort(sample(1:how.many4xbygenusgs,sample.size4xgs,replace=FALSE)) #for 6x only some 4x genera need to be sampled
	tetraploid.sample<-rep(0,sample.size4xgs)
	for(j in 1:sample.size4xgs){
		aux1<-sample(tetraploid.listgs[[sample.genus[j]]],1,replace=FALSE)
		tetraploid.sample[j]<-aux1
	}
	bootstrap.sample6x<-tetraploid.sample*1.5*rr
	bootstrapped.mean6x[i]<-mean(bootstrap.sample6x)
	bootstrapped.sd6x[i]<-sd(bootstrap.sample6x)
	
	###Creating new octoploids
	sample.genus<-sort(sample(1:how.many6xbygenusgs,sample.size6xgs,replace=FALSE)) #for 8x only some 6x genera need to be sampled
	hexaploid.sample<-rep(0,sample.size6xgs)
	for(j in 1:sample.size6xgs){
		aux1<-sample(hexaploid.listgs[[sample.genus[j]]],1,replace=FALSE)
		hexaploid.sample[j]<-aux1
	}
	
bootstrap.sample8x<-hexaploid.sample*1.33*rr
bootstrapped.mean8x[i]<-mean(bootstrap.sample8x)
bootstrapped.sd8x[i]<-sd(bootstrap.sample8x)
	
######## Creating new decaploids
	sample.genus<-sort(sample(1:how.many8xbygenusgs,sample.size8xgs,replace=FALSE)) #for 8x only some 6x genera need to be sampled
	octoploid.sample<-rep(0,sample.size8xgs)
	for(j in 1:sample.size8xgs){
		aux1<-sample(octoploid.listgs[[sample.genus[j]]],1,replace=FALSE)
		octoploid.sample[j]<-aux1
	}
	
bootstrap.sample10x<-octoploid.sample*1.25*rr
bootstrapped.mean10x[i]<-mean(bootstrap.sample10x)
bootstrapped.sd10x[i]<-sd(bootstrap.sample10x)
	

###### Creating new duodecaploids	
	decaploid.sample<-rep(0,sample.size10xgs) #for 4x all genera of 2x need to be sampled
	for(j in 1:sample.size10xgs){
		aux1<-sample(decaploid.listgs[[j]],1,replace=FALSE)
		decaploid.sample[j]<-aux1
	}

	 bootstrap.sample12x<-decaploid*1.2*rr
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

