# Descriptive statistics of chromosome number from kew.org with APGIII
# Rosana Zenil-Ferguson
#Last updated 01/27/14

# Three libraries to manipulate trees
setwd("/Users/Roszenil/Dropbox/Summer 2013/DescriptiveStats")
library(ape)
library(geiger)
library(phytools)
require(grDevices) # Color brewer
# Reading the nexus file from Soltis et al. 2011. 17 genes and mitochondrial DNA for 639 Genera. I modified it using mesquite. Original file has
# 24 trees with different genes 
#angiosperm.tree<-read.nexus(file = "familytree.nex")
# Whole tree with branch lengths. 
#plot(angiosperm.tree, edge.width=0.5,use.edge.length=FALSE, cex=0.05, no.margin=TRUE)


#Reading Dataset 
angiocvalues=read.csv("apgIIcurated2.csv")
#Defining parameters
complete<-which(angiocvalues$Chromploid==1)#5081
just.chrom<-which(angiocvalues$Minimum==1)#7121- 304 autopolyploids

# 1-Family, 2-Genus, 3-Species, 4-Authority, 5-Chromosome.number, 6-Ploidy.level, 7-Estimation.method, 8-Voucher, 9-X1C..pg., 10-Original.Reference,
# 11-Lifecycle, 12-Group, 13-Paper, 14-Name, 15-History.of.Polyploidy, 16-Repeated.Sample, 17-Non.multiple, 18-Minimum, 19-Missing.Data, 20-Possible.Hybrid, 
#21-Both, 22-Chromploid

genome.size<-angiocvalues[just.chrom,9]
fam.genome<-angiocvalues[just.chrom,1]
genus.genome<-angiocvalues[just.chrom,2]
ploidy<-angiocvalues[just.chrom,6]
chromosomes<-angiocvalues[just.chrom,5]
genera<-as.character(unique(genus.genome))
mn<-chromosomes/ploidy
mean.gs.bygenus<-tapply(genome.size,genus.genome,mean,na.rm=TRUE)
se.gs.bygenus<-tapply(genome.size, genus.genome,sd,na.rm=TRUE)
length.n.bygenus<-tapply(genome.size, genus.genome,length)

mean.mn.bygenus<-tapply(mn, genus.genome, mean, na.rm=TRUE)
sd.mn.bygenus<-tapply(mn, genus.genome,sd, na.rm=TRUE)
length.n1.bygenus<-tapply(mn, genus.genome,length)

mean.gs.byfamily<-tapply(genome.size,fam.genome,mean,na.rm=TRUE)
se.gs.byfamily<-tapply(genome.size,fam.genome,sd,na.rm=TRUE)
length.gs.byfamily<-tapply(genome.size,fam.genome,length)


sumstats.byfamily<-cbind(mean.n.byfamily, se.n.byfamily, length.n.byfamily)
write.table(sumstats.byfamily,file="mnbyfamily.cvs",sep=",",row.names=TRUE)
sumstats.byfamily<-cbind(mean.gs.byfamily,se.gs.byfamily, length.gs.byfamily)
write.table(sumstats.byfamily,file="gsbyfamily.cvs",sep=",",row.names=TRUE)

sumstatsbygenus<-cbind(mean.gs.bygenus, se.gs.bygenus,length.n.bygenus, mean.mn.bygenus,sd.mn.bygenus, length.n1.bygenus)
write.table(sumstatsbygenus,file="gsbygenus.cvs",sep=",",row.names=TRUE)
par(las=2,cex.axis=0.05)
a<-sort(length.n.bygenus[-1282])
b<-which(a==1)
length(b) #923
a<-a[-b]
barplot(height=a,width=0.01,space=1,names.arg=names(a),horiz=FALSE)

