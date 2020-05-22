#!/bin/bash

# need to get genes that are conserved in each cluster and align to get dN/dS ratios plus nt diversity

BASEDIR=/Volumes/Bennett_BACKUP/Research/jensen-metabolites/salinispora-BGCevol/BGCanalysis
BGCDIR=$BASEDIR/individualBGCs
REFDIR=$BASEDIR/evolprocesses

# going to be nearly impossible to do incomplete or split clusters
# concentrate on the better ones first

cd $REFDIR

mkdir -p $REFDIR/conservedBGCgenes
rm -rf $REFDIR/conservedBGCgenes/goodalignments
mkdir -p $REFDIR/conservedBGCgenes/goodalignments
rm -rf $REFDIR/conservedBGCgenes/badalignments
mkdir -p $REFDIR/conservedBGCgenes/badalignments

rm -f $REFDIR/conservedBGCgenes/repetitivegenes.txt
rm -f $REFDIR/conservedBGCgenes/goodconservedBGCgenes.txt

echo -e "bgcname\tcluster\tbp_count\trepseq\tgenedescription" > $REFDIR/conservedBGCgenes/goodconservedBGCgenes.txt

# for bgcname in salinosporamide
for bgcname in staurosporine sporolide salinosporamide salinipostin salinichelin lymphostin lomaiviticin salinilactam
do

	cd $BGCDIR/${bgcname}/confirmed_clusters

	rm -f *.ffn

	i=0

	for cluster in *.gbk
	do
		genomename=$(echo $cluster | cut -f1 -d'.' | sed 's/_//g')
		genbank_to_fasta.py -i $cluster -m genbank -o ${cluster%.gbk}.ffn -s nt > /dev/null 2>&1

		# add in genome name before it gets lost
		cat ${cluster%.gbk}.ffn | awk -v var="$genomename" '/^>/{sub(">", ">"var"_"++i " ")}1' > ${cluster%.gbk}.fix.ffn
		rm ${cluster%.gbk}.ffn
	done

	bgccount=$(ls *.ffn | wc -l | tr -d "[:blank:]")

	cat *.fix.ffn > total.${bgcname}.fasta
	rm -f *.fix.ffn

	cd-hit-est -i total.${bgcname}.fasta -o ${bgcname} -c 0.8 -g 1 -sc 1 > /dev/null 2>&1

	clstr2txt.pl ${bgcname}.clstr > ${bgcname}.cdhit.txt

	cat ${bgcname}.cdhit.txt | awk -v "seq=$bgccount" '{if ($3 == seq) print $0;}' > ${bgcname}.cdhit.allgenomes.txt

	cat ${bgcname}.cdhit.allgenomes.txt | cut -f2 | sort -u > ${bgcname}.cdhit.allgenomes.clusters.txt

	echo "aligning the genes found in all genomes in ${bgcname} BGC"

	while read line
	do

		cat ${bgcname}.cdhit.allgenomes.txt | awk -v "seq=$line" '{if ($2 == seq) print $0;}' | cut -f1 > ${bgcname}.cdhit.allgenomes.cluster${line}.txt
		geneIDsearch=$(cat ${bgcname}.cdhit.allgenomes.txt | awk -v "seq=$line" '{if ($2 == seq) print $0;}' | awk '{if ($5 ==1) print $1;}' )
		geneID=$(grep -w "$geneIDsearch" total.${bgcname}.fasta | cut -f2 -d'|' | awk '{$1=$1};1')
		genedescript=$(grep -w "$geneIDsearch" total.${bgcname}.fasta | cut -f3 -d'|' | awk '{$1=$1};1')

		# check to make sure each gene is only represented by one genome
		function countuniquelines () { sort "$1" | uniq -c | sort -bgr; }
		check=$(countuniquelines ${bgcname}.cdhit.allgenomes.cluster${line}.txt | head -n1 | tr -s ' ' | cut -f2 -d' ' )

		if [ $check -ne "1" ]
		then

			echo -e "${bgcname}\tcluster${line}\t${geneIDsearch}\t${geneID}" >> $REFDIR/conservedBGCgenes/repetitivegenes.txt

		else
			search-fasta.py -i total.${bgcname}.fasta -m ${bgcname}.cdhit.allgenomes.cluster${line}.txt -o total.${bgcname}.cluster${line}.temp.fasta

			sed '/^>/ s/_.*//' total.${bgcname}.cluster${line}.temp.fasta >> total.${bgcname}.cluster${line}.fasta
			clustalo -i total.${bgcname}.cluster${line}.fasta -o total.${bgcname}.cluster${line}.aln --force

			count-bp-in-seqs.sh total.${bgcname}.cluster${line}.aln > /dev/null 2>&1
			bpcount=$(grep -v ">" total.${bgcname}.cluster${line}-bp-count.txt | grep -v "sequences" | sort | uniq)

					#	# check to make sure open reading frames are good to go!
			if (( $bpcount % 3 == 0 ))           # no need for brackets
				then
					mv total.${bgcname}.cluster${line}.aln $REFDIR/conservedBGCgenes/goodalignments
				else
					mv total.${bgcname}.cluster${line}.aln $REFDIR/conservedBGCgenes/badalignments
			fi

			echo -e "${bgcname}\tcluster${line}\t${bpcount}\t${geneID}\t${genedescript}" >> $REFDIR/conservedBGCgenes/goodconservedBGCgenes.txt
			rm -f total.${bgcname}.cluster${line}.fasta
			rm -f total.${bgcname}.cluster${line}.temp.fasta
		fi

		rm -f ${bgcname}.cdhit.allgenomes.cluster${line}.txt
		rm -f total.${bgcname}.cluster${line}-bp-count.txt

	done < ${bgcname}.cdhit.allgenomes.clusters.txt

	rm -f ${bgcname}.cdhit.allgenomes.clusters.txt
	rm -f ${bgcname}.cdhit.allgenomes.txt

done

##################################################################
##################################################################
##################################################################
# need to do rif separately since it has a broken up BGC
##################################################################
##################################################################
##################################################################

##################################################################
### RIF
##################################################################

cd $BGCDIR/rifamycin/confirmed_clusters

rm -rf tempffn
mkdir -p tempffn

for f in *.gbk
do 
	genbank_to_fasta.py -i $f -m genbank -o ${f%.gbk}.ffn -s nt > /dev/null 2>&1
done

while read genome
do
	genomename=$(echo $genome | sed 's/_//g')
	cat ${genome}.*.ffn | awk -v var="$genomename" '/^>/{sub(">", ">"var"_"++i " ")}1' > tempffn/${genome}.rif.ffn
done < $BGCDIR/rifamycin/confirmed_genomes.rif.txt

rm -f *.ffn

cd $BGCDIR/rifamycin/confirmed_clusters/tempffn
bgccount=$(ls S*.ffn | cut -f1 -d'.' | wc -l | tr -d "[:blank:]")

cat *.rif.ffn > total.rif.fasta
rm -f *.ffn

cd-hit-est -i total.rif.fasta -o rif -c 0.8 -g 1 -sc 1 > /dev/null 2>&1

clstr2txt.pl rif.clstr > rif.cdhit.txt

cat rif.cdhit.txt | awk -v "seq=$bgccount" '{if ($3 == seq) print $0;}' > rif.cdhit.allgenomes.txt

cat rif.cdhit.allgenomes.txt | cut -f2 | sort -u > rif.cdhit.allgenomes.clusters.txt

echo "aligning the genes found in all genomes in rif BGC"

while read line
do

	cat rif.cdhit.allgenomes.txt | awk -v "seq=$line" '{if ($2 == seq) print $0;}' | cut -f1 > rif.cdhit.allgenomes.cluster${line}.txt
	geneIDsearch=$(cat rif.cdhit.allgenomes.txt | awk -v "seq=$line" '{if ($2 == seq) print $0;}' | awk '{if ($5 ==1) print $1;}' )
	geneID=$(grep -w "$geneIDsearch" total.rif.fasta | cut -f2 -d'|' | awk '{$1=$1};1')
	genedescript=$(grep -w "$geneIDsearch" total.rif.fasta | cut -f3 -d'|' | awk '{$1=$1};1')

	# check to make sure each gene is only represented by one genome
	function countuniquelines () { sort "$1" | cut -f1 -d'_' | uniq -c | sort -bgr; }
	check=$(countuniquelines rif.cdhit.allgenomes.cluster${line}.txt | head -n1 | tr -s ' ' | cut -f2 -d' ' )

	if [ $check -ne "1" ]
	then

		echo -e "rifamycin\tcluster${line}\t${geneIDsearch}\t${geneID}" >> $REFDIR/conservedBGCgenes/repetitivegenes.txt

	else
		search-fasta.py -i total.rif.fasta -m rif.cdhit.allgenomes.cluster${line}.txt -o total.rif.cluster${line}.temp.fasta

		sed '/^>/ s/_.*//' total.rif.cluster${line}.temp.fasta >> total.rif.cluster${line}.fasta
		clustalo -i total.rif.cluster${line}.fasta -o total.rifamycin.cluster${line}.aln --force

		count-bp-in-seqs.sh total.rifamycin.cluster${line}.aln > /dev/null 2>&1
		bpcount=$(grep -v ">" total.rifamycin.cluster${line}-bp-count.txt | grep -v "sequences" | sort | uniq)

				#	# check to make sure open reading frames are good to go!
		if (( $bpcount % 3 == 0 ))           # no need for brackets
			then
				mv total.rifamycin.cluster${line}.aln $REFDIR/conservedBGCgenes/goodalignments
			else
				mv total.rifamycin.cluster${line}.aln $REFDIR/conservedBGCgenes/badalignments
		fi

		echo -e "rifamycin\tcluster${line}\t${bpcount}\t${geneID}\t${genedescript}" >> $REFDIR/conservedBGCgenes/goodconservedBGCgenes.txt
		rm -f total.rif.cluster${line}.fasta
		rm -f total.rif.cluster${line}.temp.fasta
	fi

	rm -f rif.cdhit.allgenomes.cluster${line}.txt
	rm -f total.rifamycin.cluster${line}-bp-count.txt

done < rif.cdhit.allgenomes.clusters.txt

rm -f rif.cdhit.allgenomes.clusters.txt
rm -f rif.cdhit.allgenomes.txt

mv rif.clstr $BGCDIR/rifamycin/confirmed_clusters/rifamycin.clstr
mv rif.cdhit.txt $BGCDIR/rifamycin/confirmed_clusters/rifamycin.cdhit.txt
mv total.rif.fasta $BGCDIR/rifamycin/confirmed_clusters/total.rifamycin.fasta

rm -rf $BGCDIR/rifamycin/confirmed_clusters/tempffn


