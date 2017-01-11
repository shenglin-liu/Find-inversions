# Find inversions in the genome using high-density SNP data of multiple populations.
# Input format: VCF.
# The VCF file must be ordered according to the position of the SNPs.
# The VCF file must have 9 columns before sample columns.
# By Shenglin Liu, Nov 17, 2015.

#======================Parameters
f.vcf<-"../batch_1.vcf"	# Input file
f.popmap<-"../popmap2"	# Input file
referencePop<-"Nor"
n.digit<-4
n.nonZeroPop<-4
lowFreq<-0.1	# SNPs with overall frequency lower than this will be removed
l.window<-2e6
l.step<-1e6
fst.threshold<-0.01	# FST higher than this will be removed
minRange<-200
minInter<-200
minSnps<-3
f.output<-"Inversions.txt"

#======================Calculate allele frequency
vcf<-as.matrix(read.table(f.vcf,sep="\t",stringsAsFactors=F)[,-c(3:9)])
snpid<-paste(vcf[,1],as.integer(vcf[,2]),sep="_")
vcf<-vcf[,-c(1:2)]
nrow.vcf<-nrow(vcf)
ncol.vcf<-ncol(vcf)
chrom1<-matrix(substring(vcf,1,1),nrow.vcf,ncol.vcf)
chrom2<-matrix(substring(vcf,3,3),nrow.vcf,ncol.vcf)

popmap<-read.table(f.popmap,sep="\t",stringsAsFactors=F)[,2]
popnames<-unique(popmap)
n.pop<-length(popnames)
sampleSizes<-table(popmap)

freq<-matrix(0,nrow.vcf,n.pop)
colnames(freq)<-popnames
rownames(freq)<-snpid

for(i in 1:n.pop)
{
	index<-which(popmap==popnames[i])
	sub.chrom1<-chrom1[,index]
	sub.chrom2<-chrom2[,index]
	ref<-rowSums(sub.chrom1=="0")+rowSums(sub.chrom2=="0")
	alt<-rowSums(sub.chrom1=="1")+rowSums(sub.chrom2=="1")
	freq[,i]<-ref/(ref+alt)
}

rm(vcf,snpid,chrom1,chrom2,sub.chrom1,sub.chrom2,ref,alt)

#======================Correct allele frequency (due to missing values)
count<-freq
for(pop in popnames)	# change into count
{
	bound<-round(c(-1,seq(0,1,1/(sampleSizes[pop]*2))),n.digit)
	for(i in 0:(sampleSizes[pop]*2))
	{
		index<-((freq[,pop]>bound[i+1])&(freq[,pop]<=bound[i+2]))
		count[index,pop]<-i
	}
	freq[,pop]<-round(count[,pop]/(sampleSizes[pop]*2),n.digit)
}

rm(count)

#======================Pseudo-phase according to reference populaiton (MAF)
index<-freq[,referencePop]>0.5
freq[index,]<-1-freq[index,]
freq<-round(freq,n.digit)

#======================Filter out SNPs non-polymorphic in many populations
index<-rowSums(freq==0)<=(n.pop-n.nonZeroPop)
freq<-freq[index,]
index<-rowSums(freq==1)<=(n.pop-n.nonZeroPop)
freq<-freq[index,]

#======================Filter out SNPs with low overall allele frequency
index<-rowMeans(freq)>lowFreq
freq<-freq[index,]
index<-rowMeans(freq)<(1-lowFreq)
freq<-freq[index,]

#======================Filter out SNPs with low variance across populations
#sds<-(rowMeans(freq^2)-(rowMeans(freq))^2)^0.5
#index<-sds>mean(sds)
#freq<-freq[index,]

#======================Function for calculating FST
frq2fst<-function(freq)
{
	he<-2*freq-2*(freq^2)
	hs<-rowMeans(he)
	
	mean.freq<-rowMeans(freq)
	ht<-2*mean.freq-2*(mean.freq^2)
	
	fst<-1-hs/ht
	fst[is.na(fst)]<-0
	fst
}

#======================Function for assembling the SNPs in inverted regions
assemble<-function(possibleI)	# possibleI is a list or a dataframe
{
	assembly<-vector("list",0)
	while(T)
	{
		query<-possibleI[[1]]
		n.hit<-1
		while(T)
		{
			index<-sapply(possibleI,function(x){
				sum(query%in%x)>0
				})
			if(sum(index)==n.hit)
			{
				break
			}
			query<-unlist(possibleI[index])
			n.hit<-sum(index)
		}
		assembly<-c(assembly,list(unique(unlist(possibleI[index]))))
		possibleI<-possibleI[-which(index)]
		if(length(possibleI)==0)
		{
			break
		}
	}
	assembly
}

#======================Find inversions through sliding window
temp<-matrix(unlist(strsplit(rownames(freq),"_")),2,nrow(freq))
chr<-as.character(temp[1,])
pos<-as.integer(temp[2,])
rm(temp)
chrNames<-unique(chr)
inversions<-vector("list",0)	# A list to store the SNP positions of inversions
for(chrName in chrNames)	# Extract chromosome
{
	index<-chr==chrName
	sub.freq<-freq[index,]
	sub.pos<-pos[index]
	rownames(sub.freq)<-sub.pos
	maxPos<-max(sub.pos)
	possibleI<-character(0)
	i<-1
	while(T)	# Sliding window
	{
		ST<-(i-1)*l.step+1
		EN<-ST-1+l.window
		index<-(sub.pos>=ST)&(sub.pos<=EN)
		ss.freq<-sub.freq[index,]
		snppairs<-combn(rownames(ss.freq),2)	# SNP pairs
		fst<-numeric(0)
		for(j in 1:ncol(snppairs))	# Calculate FST for SNP pairs
		{
			fst<-cbind(fst,frq2fst(t(ss.freq[snppairs[,j],])))
		}
		fst<-colMeans(fst)
		index<-which(fst<fst.threshold)	# SNP pairs with low FST
		possibleI<-cbind(possibleI,snppairs[,index])
		if(EN>=maxPos)	# The window reaches the end of the chromosome
		{
			break
		}
		i<-i+1
	}
	sub.inversions<-assemble(data.frame(possibleI,stringsAsFactors=F))
	index1<-sapply(sub.inversions,function(x){
		(max(as.integer(x))-min(as.integer(x)))>=minRange
		})	# Filter according to range
	n.minInter<-sapply(sub.inversions,function(x){
		temp<-sort(as.integer(x))
		sum((temp[-1]-temp[-length(temp)])<=minInter)
		})	# Filter according to neighboring distance
	index2<-sapply(1:length(sub.inversions),function(x){
		(length(sub.inversions[[x]])-n.minInter[x])>=minSnps
		})	# Filter according to number of SNPs
	index<-index1&index2
	sub.inversions<-sub.inversions[index]
	inversions<-c(inversions,list(sub.inversions))
}
names(inversions)<-chrNames
index<-sapply(inversions,length)>0
inversions<-inversions[index]

#======================Output
sink(f.output)
cat(sprintf("# Input VCF: %s",f.vcf),"\n",sep="")
cat(sprintf("# Input popmap: %s",f.popmap),"\n",sep="")
cat(sprintf("# Reference population for Pseudo-phasing: %s",referencePop),"\n",sep="")
cat(sprintf("# Number of decimal digits for allele frequency: %d",n.digit),"\n",sep="")
cat(sprintf("# Least number of populations showing polymorphism for a SNP: %d",n.nonZeroPop),"\n",sep="")
cat(sprintf("# Overall allele frequency higher than: %f",lowFreq),"\n",sep="")
cat(sprintf("# Sliding window size: %d",l.window),"\n",sep="")
cat(sprintf("# Sliding step: %d",l.step),"\n",sep="")
cat(sprintf("# FST threshold: %f",fst.threshold),"\n",sep="")
cat(sprintf("# Minimum size of an inversion: %d",minRange),"\n",sep="")
cat(sprintf("# Minimum distance between neighboring SNPs: %d",minInter),"\n",sep="")
cat(sprintf("# Minimum number of SNPs for an inversion: %d",minSnps),"\n",sep="")
sink()
write("\n\n# Inversions detected:",f.output,append=T)
for(chrName in names(inversions))
{
	write(paste("\n\n# ",chrName,sep=""),f.output,append=T)
	for(i in 1:length(inversions[[chrName]]))
	{
		temp<-sort(as.integer(inversions[[chrName]][[i]]))
		write(paste("#",range(temp),collapse=" "),f.output,append=T)
		write.table(freq[paste(chrName,temp,sep="_"),],
			f.output,sep="\t",quote=F,append=T)
	}
}

#======================Filter overlap
