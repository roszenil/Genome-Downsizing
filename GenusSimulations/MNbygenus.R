# Three libraries to manipulate trees
#setwd("/Users/Roszenil/Dropbox/Summer 2013/DescriptiveStats")
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
genus.genome<-as.factor(angiocvalues[non.repeated,2]) #1631
family.genome<-as.factor(angiocvalues[non.repeated,1])#248
how.manygenus<-length(levels(genus.genome))
#minitable<-angiocvalues[non.repeated,c(2,6,9)]
#genusbyploidy=table(minitable[,1],minitable[,2])
averagesbygenus<-c(3.969318, 7.347373, 8.579877, 10.424654, 8.906167)
genus.names<-levels(genus.genome)
#dim(genusbyploidy)
#apply(genusbyploidy,2,FUN="sum")

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

	 bootstrap.sample12x<-decaploid.sample*1.2*rr
bootstrapped.mean12x[i]<-mean(bootstrap.sample12x)
bootstrapped.sd12x[i]<-sd(bootstrap.sample12x)
}

# genome.means<-#c(mean(genome.size[diploids]),mean(genome.size[tetraploids]),mean(genome.size[hexaploids]) ,mean(genome.size[octoploids]),mean(genome.size[decaploids]),mean(#genome.size[duodecaploids]))

#genome.sd<-#c(sd(genome.size[diploids]),sd(genome.size[tetraploids]),sd(genome.size[hexaploids]),sd(genome.size[octoploids]),sd(genome.size[decaploids]),sd(genome.size[d#uodecaploids]))

# Mean Squared error- Distance from observed vs bootstraped the closer to 0 the better is the prediction
qdmean4x<-sum((bootstrapped.mean4x-averagesbygenus[1])^2)/B
qdmean6x<-sum((bootstrapped.mean6x-averagesbygenus[2])^2)/B
qdmean8x<-sum((bootstrapped.mean8x-averagesbygenus[3])^2)/B
qdmean10x<-sum((bootstrapped.mean10x-averagesbygenus[4])^2)/B
qdmean12x<-sum((bootstrapped.mean12x-averagesbygenus[5])^2)/B

qdmeanh1<-c(qdmean4x,qdmean6x,qdmean8x,qdmean10x,qdmean12x) # Vector of distances
return(qdmeanh1)
}

B<-10000 #number of bootstraps
qdmean<-retention.rate2(rr=1,B)

write.table(qdmean,file="bootstrapbygenusmn.txt",sep="\n",row.names=FALSE,col.names=FALSE)

## MN+GD
ret.rates<-seq(0.001,1,0.001) # retention rates tested
long1<-length(ret.rates)
 #Result of each hypothesis
bootstrapped.sim1<-rep(0,5)# Samples with each rate
for(j in 1:long1){
	qdmean<-retention.rate2(rr=ret.rates[j],B)
	bootstrapped.sim1<-rbind(bootstrapped.sim1,qdmean)	
}

write.table(bootstrapped.sim1,file="bootstrapbygenusmngd.txt",sep="\n",row.names=FALSE,col.names=FALSE)



