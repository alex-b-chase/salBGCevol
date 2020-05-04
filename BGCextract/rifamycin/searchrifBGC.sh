#!/bin/bash

# look for lymA and lymB genes for the lymphostin BGC
# OR
# look for salA and salB genes for the salinosporamide BGC
# OR
# look for kinA and kinB genes for the lomaiviticin BGC
# OR
# look for kinA and kinB genes for the staurosporine BGC
# OR
# look for kinA and kinB genes for the salinipostin BGC

gene=rif
genedir=rifamycin
genefile=${gene}EFgenes.faa

BASEDIR=/Volumes/Bennett_BACKUP/Research/jensen-metabolites/salinispora-BGCevol

GENDIR=$BASEDIR/genomic-data/prokka_annotation
REFDIR=$BASEDIR/bigscape
OUTDIR=$BASEDIR/BGCanalysis/individualBGCs/$genedir

# going to use the BIGSCAPE network to get the antismash clusters and parse those
# what we want to do is use the core biosynthesis genes 
# run this in R and can get the cluster involved in network, then screen for core biosynthesis genes

cd $OUTDIR

rm -rf $OUTDIR/confirmed_clusters
mkdir -p $OUTDIR/confirmed_clusters

grep -v "BGC" ${gene}.module.txt > ${gene}.module.temp.txt

while read line
do

	cp $REFDIR/bigscape_out/cache/fasta/${line}.fasta ./${line}.faa

done < ${gene}.module.temp.txt

rm -f *.gbk
rm -f ${gene}.module.temp.txt
rm -f $OUTDIR/${gene}EF.percid.txt
rm -f confirmed_clusters.${gene}.txt

# number of genes screening for
check=2
# distance between genes (i.e. 1 = sequential genes)
check2=1

for f in S*.faa
do 
	genome=$(echo $f | cut -f1 -d'.')
	contig=$(echo $f | cut -f2 -d'.' | sed "s/contig_/c/g")
	region=$(echo $f | cut -f3 -d'.' | sed "s/region/-/g")
	cluster=${contig}${region}
	filebasename=${f%.faa}
	echo "$genome"
	blat -prot -out=blast8 $genefile $f ${genome}.${cluster}.blat.txt

	cat ${genome}.${cluster}.blat.txt | sort -k2,2 -k4,4nr -k3,3nr | sort -u -k2,2 > ${genome}.${cluster}.blat.top.txt
	count=$(cat ${genome}.${cluster}.blat.top.txt | wc -l)
	echo ""

	if [ $count -eq $check ]
	then

		# now make sure they are sequential
		num1=$(cut -f1 -d':' ${genome}.${cluster}.blat.top.txt | rev | cut -f1 -d'_' | rev | sed 's/ORF//g' | sort -un | head -n1)
		num2=$(cut -f1 -d':' ${genome}.${cluster}.blat.top.txt | rev | cut -f1 -d'_' | rev | sed 's/ORF//g' | sort -un | tail -n1)
		final=$(expr $num2 - $num1)

		if [ $final -eq $check2 ]
		then

			cp $REFDIR/input_files/${filebasename}* $OUTDIR/confirmed_clusters/

			echo $filebasename >> confirmed_clusters.${gene}.txt

			cat ${genome}.${cluster}.blat.top.txt | grep "${gene}E" | cut -f1 > $genome.tempE.txt
			cat ${genome}.${cluster}.blat.top.txt | grep "${gene}F" | cut -f1 > $genome.tempF.txt
			search-fasta.py -i $f -m $genome.tempE.txt -o $genome.${gene}E.temp.fasta
			search-fasta.py -i $f -m $genome.tempF.txt -o $genome.${gene}F.temp.fasta

			cat $genome.${gene}E.temp.fasta | sed "s/[>].*$/>$genome/" > $genome.${gene}E.fasta
			cat $genome.${gene}F.temp.fasta | sed "s/[>].*$/>$genome/" > $genome.${gene}F.fasta

			percidA=$(cat ${genome}.${cluster}.blat.top.txt | grep "${gene}E" | cut -f3 )
			percidB=$(cat ${genome}.${cluster}.blat.top.txt | grep "${gene}F" | cut -f3 )

			echo -e "${genome}\t${gene}E\t${percidA}" >> $OUTDIR/${gene}EF.percid.txt
			echo -e "${genome}\t${gene}F\t${percidB}" >> $OUTDIR/${gene}EF.percid.txt


			rm -f $genome.tempE.txt
			rm -f $genome.tempF.txt
			rm -f $genome.${gene}E.temp.fasta
			rm -f $genome.${gene}F.temp.fasta
		else
			:
		fi
	else
		:
	fi

	rm -f ${genome}.${cluster}.blat.txt
	rm -f ${genome}.${cluster}.blat.top.txt
	rm -f $f 

done


cat *.${gene}E.fasta > $OUTDIR/${gene}E.total.fasta
cat *.${gene}F.fasta > $OUTDIR/${gene}F.total.fasta

rm -f *.${gene}E.fasta 
rm -f *.${gene}F.fasta 

cd $OUTDIR

clustalo -i ${gene}E.total.fasta -o ${gene}E.total.aln --force -v
clustalo -i ${gene}F.total.fasta -o ${gene}F.total.aln --force -v

rm -f ${gene}E.total.fasta
rm -f ${gene}F.total.fasta

catfasta2phyml.pl -f *.aln  > ${gene}.aligned.filtered.fa

rm -f RAxML_*.${gene}

raxml -s ${gene}.aligned.filtered.fa -m PROTGAMMAWAG -n ${gene} -x 100 -# 100 -p 4321 -f a -T 4

mkdir -p $OUTDIR/$gene.EFphylo
mv *.aln $OUTDIR/$gene.EFphylo
mv ${gene}.aligned.filtered.fa* $OUTDIR/$gene.EFphylo
mv RAxML_* $OUTDIR/$gene.EFphylo


