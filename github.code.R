#In this script, I include code for demultiplexing, trimming, generating multi-locus genotypes, and running population genomic analyses
#This code includes examples of each step in the Linux bioinformatic and R pipeline 
#This code was produced to study the population genomics of a morphologically convergent coral species group called massive Porites
#Found throughout the Indo-Pacific, they occur in different reef environments, and in this study, we find that different species,
#not populations of the same species, predominantly occupy river deltas and adjacent (< 100M) fore reefs on Guam respectively, 
#and bleaching less prevalent in river deltas. 

#Demultiplexing fastq files in order to generate sample-specific fastq files
screen -L python3 identify_dbrs6.py \n
-i ~/Karim_Porites/ddRAD/info/sampleinfofiles/RAD.b0.sampleinfofiles/RAD8.3.updated.sampleinfofile.txt -l7 1 \n
-f ~/Karim_Porites/ddRAD/raw/RAD8.3.trimming.output -s KP.RAD8.3.demult.output.txt -b 0

#Building bowtie reference
bowtie2-build /media/RAID/david/Genome/Porites_lutea/BWA-host/Plut.fasta.gz Porites_lutea.fasta.gz
bowtie2-build /media/RAID/david/Genome/Porites_lutea/BWA-symb/SymbC15_plutea_v2.1.fna symbC15.plutea.fasta.gz
bowtie2-build /media/RAID/david/Genome/Porites_lutea/BWA-bact/plut_bact.fa plut_bact.fasta.gz


#Using bowtie to align reads to reference genomes (first the symbionts and microbiome to remove microbe and symbiont reads,
#then aligning to the host genome to recover host reads)

#!/bin/bash
for sample in ./*.fq.gz
do
base=$(basename "$sample" .fq.gz)
bowtie2 -x ~/Karim_Porites/ddRAD/bowtie/bowtie.coral.reference/P.mass.coral.reference -U ${base}.fq.gz --un ~/Karim_Porites/ddRAD/bowtie/coral.nosplit.techreps/coral.un/${base}_trash.sam --al ~/Karim_Porites/ddRAD/bowtie/coral.nosplit.techreps/coral.al/${base}_keep.fq.gz -S ~/Karim_Porites/ddRAD/bowtie/coral.nosplit.techreps/coral.al/${base}_keep.sam 2>>bowtie.coralfqtosam.nosplit.techreps.log
done

#Symbiont (keep un, trash al!!!)
Pwd for input data: ~/Karim_Porites/ddRAD/bowtie/bacteria.nosplit.techreps/bacteria.un
#!/bin/bash
for sample in ./*.fq.gz
do
base=$(basename "$sample" .fq.gz)
bowtie2 -x ~/Karim_Porites/ddRAD/bowtie/bowtie.symbiont.reference/P.mass.symbiont.reference -U ${base}.fq.gz --un ~/Karim_Porites/ddRAD/bowtie/symbiont.nosplit.techreps/symbiont.un/${base}_keep.fq.gz --al ~/Karim_Porites/ddRAD/bowtie/symbiont.nosplit.techreps/symbiont.al/${base}_trash.fq.gz -S ~/Karim_Porites/ddRAD/bowtie/symbiont.nosplit.techreps/symbiont.al/${base}_trash.sam 2>>bowtie.symbiont.nosplit.techreps.log
done

#Bacteria (keep un, trash al!!!)
Pwd for input data: ~/Karim_Porites/ddRAD/cleaned/thesis.nosplit.techreps.200806.dataset.gzipped
#!/bin/bash
for sample in ./*.fq.gz
do
base=$(basename "$sample" .fq.gz)
bowtie2 -x ~/Karim_Porites/ddRAD/bowtie/bowtie.bacteria.reference/P.mass.bacteria.reference -U ${base}.fq.gz --un ~/Karim_Porites/ddRAD/bowtie/bacteria.nosplit.techreps/bacteria.un/${base}_keep.fq.gz --al ~/Karim_Porites/ddRAD/bowtie/bacteria.nosplit.techreps/bacteria.al/${base}_trash.fq.gz 2>>bowtie.bacteria.nosplit.techreps.log
Done

#convert sam to bam files
#Sam → bam loop with -F 260 flag
#!/bin/bash
for i in ./*.coral.sam
do 
samtools view -S -b -F 260 "$i" > "${i%.coral.sam}".bam 
done

#Sort bam files 
#!/bin/bash
for i in ./*.bam
do 
samtools sort "$i" > "${i%.bam}".sorted.bam 
done

#Index bam files 
for i in *.sorted.bam
do 
samtools index $i > ${i%.sorted.bam}.indexed.bam
Done


#Phylogenetic Tree Code
#Phylogenetic tree phylip file code:
/usr/local/bin/populations -P . -M /media/RAID/karim/Karim_Porites/ddRAD/info/popmaps/thesis.dataset.popmaps/200815.bowtie.aligned.phylogenetic.popmap -R 0.30 --phylip-var


#RAxML Phylogenetic Tree Reconstruction Code:
/usr/bin/raxmlHPC-PTHREADS-AVX -T 4 -n result -s 200829.10bs.iqtree.varsites.phy -m ASC_GTRCAT -c 25 -p 12345 -f a -N 1000 -x 12345 --asc-corr lewis


Fastq to Genotype Likelihood Analysis

#Examples:

#To generate ibs, cov.Mat, and other GL files:
#Genomic dataset:
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 79 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -setMinDepthInd 8 -doCov 1 -doGeno 8 -doVcf 1 -doPost 1 -doGlf 2"
/media/RAID/david/tools/ANGSD/angsd/angsd -b wholedataset.heterozygosity.bams -GL 1 $FILTERS $TODO -P 1 -out 201014.wholedataset.heterozygosity
Pink clade:
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 10 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -setMinDepthInd 8 -doCov 1 -doGeno 8 -doVcf 1 -doPost 1 -doGlf 2"
/media/RAID/david/tools/ANGSD/angsd/angsd -b pinkclade.pop.bams -GL 1 $FILTERS $TODO -P 1 -out pinkclade.ibs05

#To generate beagle files (which are used as input for NGSAdmix plots):
#Orange clade:
/media/RAID/david/tools/ANGSD/angsd/angsd -GL 1 -out genolike -nThreads 10 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -bam orangeclade.pop.bams


#PCoA Example
#ANGSD:
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 16 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -setMinDepthInd 8 -doCov 1 -doGeno 8 -dobcf 1 -doPost 1 -doGlf 2"
/media/SCRATCH/hector/0SOFTWARE/ANGSD/angsd -b 230303.clades1.2.4.genomic.dataset.admix.bams -GL 1 $FILTERS $TODO -P 1 -out 230303.clades1.2.4.genomic.datasetpcoa.ibs05

scp ~/Desktop/210302.publication.genomicdatasetpcoa.bams/ karim@168.123.185.34:/media/RAID/karim/Karim_Porites/ddRAD/bowtie.angsd.nosplit/bowtie.nonsplit.final.popgen.btwnspecies.dataset/publication.dataset/ 
scp karim@168.123.185.34:/media/RAID/karim/Karim_Porites/ddRAD/bowtie.angsd.nosplit/bowtie.nonsplit.final.popgen.btwnspecies.dataset/publication.dataset/210302.publication.genomicdatasetpcoa.i2p ~/Desktop/
scp karim@168.123.185.34:/media/RAID/karim/Karim_Porites/ddRAD/bowtie.angsd.nosplit/bowtie.nonsplit.final.popgen.btwnspecies.dataset/publication.dataset/210302.publication.genomicdatasetpcoa.ibs05.covMat ~/Desktop/
  
  
#In R Studio  
library(Rtsne)
library(vegan)
library(adegenet)
quartz ()
bams=read.table("~/Desktop/210302.publication.genomicdatasetpcoa.bams")[,1]
goods=c(1:length(bams))
i2p=read.table("~/Desktop/210302.publication.genomicdatasetpcoa.i2p", sep="\t")
row.names(i2p)=i2p[,1]
i2p=i2p[goods,]
site=i2p[,2]
enviro=i2p[,3]
clade=i2p[,4]
shapes=c(16,17)
shapes <-shapes[as.numeric(as.factor(enviro))]
palette(c("blue","green", "orange", "red"))
colors=as.numeric(as.factor(site))
colpops=as.numeric(as.factor(sort(unique(site))))
co=as.matrix(read.table("~/Desktop/210302.publication.genomicdatasetpcoa.ibs05.covMat"))
co=co[goods,goods]
dimnames(co)=list(bams[goods],bams[goods])
conds=data.frame(cbind(site))
pp0=capscale(as.dist(1-cov2cor(co))~1)
axes2plot=c(1,2)
quartz()
cc=pp0
plot(cc, choices=axes2plot,type="n")
points(cc, choices=axes2plot, pch=shapes, col=colors)
ordispider(cc, choices=axes2plot,groups=clade,col="grey80", label=T)
legend("topleft", legend = c("Fouha","Inarajan", "Ritidian", "Talofofo"), col = colpops, pch = 19, bty = "7", pt.cex = 2, cex = 1.2, text.col = "black", horiz = F , inset = c(0.1, 0.1))
legend("bottomleft", legend=c("Fore Reef", "River Delta"), col= "black", pch = c(16, 17), bty = "7", pt.cex = 2, cex = 1.2, text.col = "black", horiz = F , inset = c(0.1, 0.1))



#Admixture Plot Example
(K from 2 - 8)
for K in `seq 2 8`
do
~/tools/NGSadmix -likes 230303.clades4to7.genomic.dataset.admixture.beagle.gz -K $K -P 10 -o 230303.clades4to7.genomic.dataset.ngsadmix_k${K}
done

scp karim@168.123.185.34:/media/RAID/karim/Karim_Porites/ddRAD/bowtie.angsd.nosplit/bowtie.nonsplit.final.popgen.btwnspecies.dataset/publication.dataset/210302.publication.genomicdataset.admix.bams ~/Desktop/
#Run on cluster, not R!
R --vanilla --slave -e 'admix <- pop<-read.table("clades.4to7.pop.info",as.is=T);q<-read.table("230303.clades4to7.genomic.dataset.ngsadmix_k4.qopt");ord<-order(pop[,1]);barplot(t(q)[,ord],col=1:10,space=0.01,border=NA,xlab="Clade",ylab="Admixture Proportions for K=4");text(tapply(1:nrow(pop),pop[ord,1],mean),-0.05,unique(pop[ord,1]),xpd=T);abline(v=cumsum(sapply(unique(pop[ord,1]),function(x){sum(pop[ord,1]==x)})),col=1,lwd=1.2)'

scp karim@168.123.185.34:/media/RAID/karim/Karim_Porites/ddRAD/bowtie.angsd.nosplit/bowtie.nonsplit.final.popgen.btwnspecies.dataset/publication.dataset/230303.clades4to7.genomic.dataset.admix/clades.4to7.pop.info .
scp karimp@shell.uog.edu:/home/karimp/clades.4to7.pop.info ~/Desktop/
scp karim@168.123.185.34:/media/RAID/karim/Karim_Porites/ddRAD/bowtie.angsd.nosplit/bowtie.nonsplit.final.popgen.btwnspecies.dataset/publication.dataset/230303.clades4to7.genomic.dataset.admix/230303.clades4to7.genomic.dataset.ngsadmix_k2.qopt .
scp karimp@shell.uog.edu:/home/karimp/230303.clades4to7.genomic.dataset.ngsadmix_k4.qopt ~/Desktop/
  
  










































