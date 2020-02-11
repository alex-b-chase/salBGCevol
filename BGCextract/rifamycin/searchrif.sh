#!/bin/bash

# look for lymA and lymB genes for the lymphostin BGC
# OR
# look for salA and salB genes for the salinosporamide BGC

gene=rif
genedir=rifamycin_search
genefile=${gene}APgenes.faa
genefile2=${gene}EFgenes.faa

BASEDIR=/Volumes/Bennett_BACKUP/Research/jensen-metabolites/salinispora-genomics
GENDIR=$BASEDIR/prodigal_annotation
REFDIR=$BASEDIR/antismash/bigscape
OUTDIR=$REFDIR/$genedir

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

# genome bad!
rm -f SA_CNY281*
rm -f $OUTDIR/${gene}DE.percid.txt
rm -f $OUTDIR/${gene}AP.percid.txt
rm -f confirmed_clusters.${gene}.txt

check=2
check2=1

for f in S*.faa
do 
	genome=$(echo $f | cut -f1 -d'.' | sed 's/_//g')
	cluster=${f%.faa}
	echo "$genome"
	blat -prot -out=blast8 $genefile2 $f ${cluster}.blat.txt

	cat ${cluster}.blat.txt | sort -k2,2 -k4,4nr -k3,3nr | sort -u -k2,2 > ${genome}.blat.top.txt
	count=$(cut -f1 -d':' ${genome}.blat.top.txt | sort -u | wc -l)
	echo ""

	if [ $count -eq $check ]
	then

		# now make sure they are sequential
		num1=$(cut -f1 -d':' ${genome}.blat.top.txt | cut -f3 -d'_' | sed 's/ORF//g' | sort -un | head -n1)
		num2=$(cut -f1 -d':' ${genome}.blat.top.txt | cut -f3 -d'_' | sed 's/ORF//g' | sort -un | tail -n1)
		final=$(expr $num2 - $num1)

		if [ $final -eq $check2 ]
		then

			cp $REFDIR/input_files/${cluster}* $OUTDIR/confirmed_clusters/

			echo $cluster >> confirmed_clusters.${gene}.txt

			cat ${genome}.blat.top.txt | grep "${gene}F" | cut -f1 > $genome.tempF.txt
			cat ${genome}.blat.top.txt | grep "${gene}E" | cut -f1 > $genome.tempE.txt
			search-fasta.py -i $f -m $genome.tempF.txt -o $genome.${gene}F.temp.fasta
			search-fasta.py -i $f -m $genome.tempE.txt -o $genome.${gene}E.temp.fasta

			cat $genome.${gene}F.temp.fasta | sed "s/[>].*$/>$genome/" > $genome.${gene}F.fasta
			cat $genome.${gene}E.temp.fasta | sed "s/[>].*$/>$genome/" > $genome.${gene}E.fasta

			percidD=$(cat ${genome}.blat.top.txt | grep "${gene}F" | cut -f3 )
			percidE=$(cat ${genome}.blat.top.txt | grep "${gene}E" | cut -f3 )

			echo -e "${genome}\t${gene}F\t${percidD}" >> $OUTDIR/${gene}EF.percid.txt
			echo -e "${genome}\t${gene}E\t${percidE}" >> $OUTDIR/${gene}EF.percid.txt


			rm -f $genome.tempF.txt
			rm -f $genome.tempE.txt
			rm -f $genome.${gene}F.temp.fasta
			rm -f $genome.${gene}E.temp.fasta
		else
			:
		fi
	elif [ $count -lt $check ]
	then
		# check if the genes are the rifA and P450 genes for rifBGC	
		blat -prot -out=blast8 $genefile $f ${cluster}.blat.txt

		cat ${cluster}.blat.txt | sort -k2,2 -k4,4nr -k3,3nr | sort -u -k2,2 > ${genome}.blat.top.txt
		count2=$(cut -f1 -d':' ${genome}.blat.top.txt | sort -u | wc -l)

		if [ $count2 -eq $check ]
		then

			# now make sure they are sequential
			num1=$(cut -f1 -d':' ${genome}.blat.top.txt | cut -f3 -d'_' | sed 's/ORF//g' | sort -un | head -n1)
			num2=$(cut -f1 -d':' ${genome}.blat.top.txt | cut -f3 -d'_' | sed 's/ORF//g' | sort -un | tail -n1)
			final=$(expr $num2 - $num1)

			if [ $final -lt 3 ]
			then

				cp $REFDIR/input_files/${cluster}* $OUTDIR/confirmed_clusters/

				echo $cluster >> confirmed_clusters.${gene}.txt

				cat ${genome}.blat.top.txt | grep "${gene}A" | cut -f1 > $genome.tempA.txt
				cat ${genome}.blat.top.txt | grep "${gene}P" | cut -f1 > $genome.tempP.txt
				search-fasta.py -i $f -m $genome.tempA.txt -o $genome.${gene}A.temp.fasta
				search-fasta.py -i $f -m $genome.tempP.txt -o $genome.${gene}P.temp.fasta

				cat $genome.${gene}A.temp.fasta | sed "s/[>].*$/>$genome/" > $genome.${gene}A.fasta
				cat $genome.${gene}P.temp.fasta | sed "s/[>].*$/>$genome/" > $genome.${gene}P.fasta

				percidA=$(cat ${genome}.blat.top.txt | grep "${gene}A" | cut -f3 )
				percidP=$(cat ${genome}.blat.top.txt | grep "${gene}P" | cut -f3 )

				echo -e "${genome}\t${gene}A\t${percidA}" >> $OUTDIR/${gene}AP.percid.txt
				echo -e "${genome}\t${gene}P\t${percidP}" >> $OUTDIR/${gene}AP.percid.txt


				rm -f $genome.tempA.txt
				rm -f $genome.tempP.txt
				rm -f $genome.${gene}A.temp.fasta
				rm -f $genome.${gene}P.temp.fasta
			else
				:
			fi
		else
			:
		fi
	else
		:
	fi

	rm -f ${cluster}.blat.txt
	rm -f ${genome}.blat.top.txt
	rm -f $f 

done


cat *.${gene}F.fasta > $OUTDIR/${gene}F.total.fasta
cat *.${gene}E.fasta > $OUTDIR/${gene}E.total.fasta

cat *.${gene}A.fasta > $OUTDIR/${gene}A.total.fasta
cat *.${gene}P.fasta > $OUTDIR/${gene}P.total.fasta

rm -f *.${gene}F.fasta 
rm -f *.${gene}E.fasta 
rm -f *.${gene}A.fasta 
rm -f *.${gene}P.fasta 

cd $OUTDIR

clustalo -i ${gene}F.total.fasta -o ${gene}F.total.aln --force -v
clustalo -i ${gene}E.total.fasta -o ${gene}E.total.aln --force -v

catfasta2phyml.pl -f *.aln  > ${gene}.aligned.filtered.fa

# rm -f RAxML_*.${gene}

# raxml -s ${gene}.aligned.filtered.fa -m PROTGAMMAWAG -n ${gene} -x 100 -# 100 -p 4321 -f a -T 4



# cd $REFDIR

# source activate bigscape

# BIGDIR=/anaconda2/envs/bigscape/BiG-SCAPE
# PFAMDIR=/Users/alexchase/software/PfamScan

# python $BIGDIR/bigscape.py -i input_files -o rifamycin \
# --pfam_dir $PFAMDIR -c 6 --mibig --cutoffs 0.3 \
# --query_bgc input_files/BGC0000137.1.fasta --clan_cutoff 0.3 0.9 \
# --mix --no_classify --mode global

# source deactivate bigscape
