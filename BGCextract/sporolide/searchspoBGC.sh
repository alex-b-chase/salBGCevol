#!/bin/bash

gene=spo
genedir=sporolide
genefile=${gene}EPgenes.faa

BASEDIR=/Volumes/Bennett_BACKUP/Research/jensen-metabolites/salinispora-BGCevol

GENDIR=$BASEDIR/genomic-data/prokka_annotation
REFDIR=$BASEDIR/bigscape
OUTDIR=$BASEDIR/BGCanalysis/individualBGCs/$genedir

################################ OLD STUFF ################################
# cd $GENDIR

# find . -name "*.ffn" -exec cp {} ../antismash/bigscape/lymphostin_search/ \;

# run BLAT againt reference BGC but need to account for alignment length
# going to subset the best match based on alignment length and then get best percID
################################ OLD STUFF ################################

# going to use the BIGSCAPE network to get the antismash clusters and parse those
# what we want to do is use the core biosynthesis genes 
# run this in R and can get the cluster involved in network, then screen for core biosynthesis genes

cd $OUTDIR

mkdir -p $OUTDIR/confirmed_clusters

grep -v "BGC" ${gene}.module.txt >${gene}.module.temp.txt

while read line
do

	cp $REFDIR/bigscape_out/cache/fasta/${line}.fasta ./${line}.faa

done < ${gene}.module.temp.txt

rm -f *.gbk
rm -f ${gene}.module.temp.txt
rm -f $OUTDIR/${gene}AB.percid.txt
rm -f confirmed_clusters.${gene}.txt

check=2
check2=3

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

			cat ${genome}.${cluster}.blat.top.txt | grep "${gene}E" | cut -f1 > $genome.tempA.txt
			cat ${genome}.${cluster}.blat.top.txt | grep "${gene}P" | cut -f1 > $genome.tempB.txt
			search-fasta.py -i $f -m $genome.tempA.txt -o $genome.${gene}A.temp.fasta
			search-fasta.py -i $f -m $genome.tempB.txt -o $genome.${gene}B.temp.fasta

			cat $genome.${gene}A.temp.fasta | sed "s/[>].*$/>$genome/" > $genome.${gene}A.fasta
			cat $genome.${gene}B.temp.fasta | sed "s/[>].*$/>$genome/" > $genome.${gene}B.fasta

			percidA=$(cat ${genome}.${cluster}.blat.top.txt | grep "${gene}E" | cut -f3 )
			percidB=$(cat ${genome}.${cluster}.blat.top.txt | grep "${gene}P" | cut -f3 )

			echo -e "${genome}\t${gene}E\t${percidA}" >> $OUTDIR/${gene}EP.percid.txt
			echo -e "${genome}\t${gene}P\t${percidB}" >> $OUTDIR/${gene}EP.percid.txt


			rm -f $genome.tempA.txt
			rm -f $genome.tempB.txt
			rm -f $genome.${gene}A.temp.fasta
			rm -f $genome.${gene}B.temp.fasta
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


cat *.${gene}A.fasta > $OUTDIR/${gene}A.total.fasta
cat *.${gene}B.fasta > $OUTDIR/${gene}B.total.fasta

rm -f *.${gene}A.fasta 
rm -f *.${gene}B.fasta 

cd $OUTDIR

clustalo -i ${gene}A.total.fasta -o ${gene}E.total.aln --force -v
clustalo -i ${gene}B.total.fasta -o ${gene}P.total.aln --force -v

rm -f ${gene}A.total.fasta
rm -f ${gene}B.total.fasta

catfasta2phyml.pl -f *.aln  > ${gene}.aligned.filtered.fa

rm -f RAxML_*.${gene}

raxml -s ${gene}.aligned.filtered.fa -m PROTGAMMAWAG -n ${gene} -x 100 -# 100 -p 4321 -f a -T 4

mkdir -p $OUTDIR/$gene.EPphylo
mv *.aln $OUTDIR/$gene.EPphylo
mv ${gene}.aligned.filtered.fa* $OUTDIR/$gene.EPphylo
mv RAxML_* $OUTDIR/$gene.EPphylo


# cd $REFDIR

# source activate bigscape

# BIGDIR=/anaconda2/envs/bigscape/BiG-SCAPE
# PFAMDIR=/Users/alexchase/software/PfamScan

# python $BIGDIR/bigscape.py -i input_files -o rifamycin \
# --pfam_dir $PFAMDIR -c 6 --mibig --cutoffs 0.3 \
# --query_bgc input_files/BGC0000137.1.fasta --clan_cutoff 0.3 0.9 \
# --mix --no_classify --mode global

# source deactivate bigscape
