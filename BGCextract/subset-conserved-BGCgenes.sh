#!/bin/bash

# need to get genes that are conserved in each cluster and align to get dN/dS ratios plus nt diversity

REFDIR=/Volumes/Bennett_BACKUP/Research/jensen-metabolites/salinispora-genomics/antismash/bigscape

# going to be nearly impossible to do incomplete or split clusters
# concentrate on the better ones first

cd $REFDIR

# rm -rf $REFDIR/coreBGCgenes
mkdir -p $REFDIR/coreBGCgenes
rm -rf $REFDIR/coreBGCgenes/goodalignments
mkdir -p $REFDIR/coreBGCgenes/goodalignments
rm -rf $REFDIR/coreBGCgenes/badalignments
mkdir -p $REFDIR/coreBGCgenes/badalignments

rm -f $REFDIR/coreBGCgenes/repetitivegenes.txt
rm -f $REFDIR/coreBGCgenes/goodcoreBGCgenes.txt

echo -e "bgcname\tcluster\tbp_count\trepseq\tgenedescription" > $REFDIR/coreBGCgenes/goodcoreBGCgenes.txt

# for bgcname in salinosporamide
for bgcname in staurosporine sporolide salinosporamide salinipostin salinichelin lymphostin lomaiviticin cyanosporaside
do

	cd $REFDIR/${bgcname}_search/confirmed_clusters

	rm -f *.ffn
	rm -f total.*.cluster*.fasta

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

			echo -e "${bgcname}\tcluster${line}\t${geneIDsearch}\t${geneID}" >> $REFDIR/coreBGCgenes/repetitivegenes.txt

		else
			search-fasta.py -i total.${bgcname}.fasta -m ${bgcname}.cdhit.allgenomes.cluster${line}.txt -o total.${bgcname}.cluster${line}.temp.fasta

			sed '/^>/ s/_.*//' total.${bgcname}.cluster${line}.temp.fasta >> total.${bgcname}.cluster${line}.fasta
			clustalo -i total.${bgcname}.cluster${line}.fasta -o total.${bgcname}.cluster${line}.aln --force

			count-bp-in-seqs.sh total.${bgcname}.cluster${line}.aln > /dev/null 2>&1
			bpcount=$(grep -v ">" total.${bgcname}.cluster${line}-bp-count.txt | grep -v "sequences" | sort | uniq)

					#	# check to make sure open reading frames are good to go!
			if (( $bpcount % 3 == 0 ))           # no need for brackets
				then
					mv total.${bgcname}.cluster${line}.aln $REFDIR/coreBGCgenes/goodalignments
				else
					mv total.${bgcname}.cluster${line}.aln $REFDIR/coreBGCgenes/badalignments
			fi

			echo -e "${bgcname}\tcluster${line}\t${bpcount}\t${geneID}\t${genedescript}" >> $REFDIR/coreBGCgenes/goodcoreBGCgenes.txt
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
# need to do rif and slm separately since they have broken up BGCs
##################################################################
##################################################################
##################################################################

##################################################################
### RIF
##################################################################

cd $REFDIR/rifamycin_search

rm -rf tempffn
mkdir -p tempffn
cd $REFDIR/rifamycin_search/tempffn/

cat $REFDIR/rifamycin_search/rif.module.txt | grep -v "BGC" > rif.module.txt

# remove incomplete rif BGCs


while read cluster
do

	genomename=$(echo $cluster | cut -f1 -d'.' | sed 's/_//g')

	if grep -Fxq "$genomename" $REFDIR/rifamycin_search/confirmed_clusters/incomplete.BGCs.txt
	then
		:
	else
		cp $REFDIR/input_files/${cluster}* $REFDIR/rifamycin_search/tempffn/
		genbank_to_fasta.py -i ${cluster}.gbk -m genbank -o ${cluster}.ffn -s nt > /dev/null 2>&1
	fi

done < rif.module.txt

# combine files from same genome
ls S*.ffn | cut -f1 -d'.' | sort -u > genomes.txt

while read genome 
do
	genomename=$(echo $genome | sed 's/_//g')
	# add in genome name before it gets lost
	cat ${genome}.cluster*.ffn | awk -v var="$genomename" '/^>/{sub(">", ">"var"_"++i " ")}1' > ${genome}.fix.ffn
	rm ${genome}.cluster*.ffn
done < genomes.txt

# now start the clustering step
cd $REFDIR/rifamycin_search/confirmed_clusters
bgccount=$(ls S*.fasta | cut -f1 -d'.' | wc -l | tr -d "[:blank:]")

cd $REFDIR/rifamycin_search/tempffn/
cat *.fix.ffn > total.rif.fasta
rm -f *.fix.ffn
rm -f *.gbk

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

		echo -e "rifamycin\tcluster${line}\t${geneIDsearch}\t${geneID}" >> $REFDIR/coreBGCgenes/repetitivegenes.txt

	else
		search-fasta.py -i total.rif.fasta -m rif.cdhit.allgenomes.cluster${line}.txt -o total.rif.cluster${line}.temp.fasta

		sed '/^>/ s/_.*//' total.rif.cluster${line}.temp.fasta >> total.rif.cluster${line}.fasta
		clustalo -i total.rif.cluster${line}.fasta -o total.rifamycin.cluster${line}.aln --force

		count-bp-in-seqs.sh total.rifamycin.cluster${line}.aln > /dev/null 2>&1
		bpcount=$(grep -v ">" total.rifamycin.cluster${line}-bp-count.txt | grep -v "sequences" | sort | uniq)

				#	# check to make sure open reading frames are good to go!
		if (( $bpcount % 3 == 0 ))           # no need for brackets
			then
				mv total.rifamycin.cluster${line}.aln $REFDIR/coreBGCgenes/goodalignments
			else
				mv total.rifamycin.cluster${line}.aln $REFDIR/coreBGCgenes/badalignments
		fi

		echo -e "rifamycin\tcluster${line}\t${bpcount}\t${geneID}\t${genedescript}" >> $REFDIR/coreBGCgenes/goodcoreBGCgenes.txt
		rm -f total.rif.cluster${line}.fasta
		rm -f total.rif.cluster${line}.temp.fasta
	fi

	rm -f rif.cdhit.allgenomes.cluster${line}.txt
	rm -f total.rifamycin.cluster${line}-bp-count.txt

done < rif.cdhit.allgenomes.clusters.txt

rm -f rif.cdhit.allgenomes.clusters.txt
rm -f rif.cdhit.allgenomes.txt



##################################################################
### SLM
##################################################################
cd $REFDIR/salinilactam_search

rm -rf tempffn
mkdir -p tempffn
cd $REFDIR/salinilactam_search/tempffn/

while read cluster
do
	genome=$(echo $cluster | cut -f1 -d'.')

	if grep -Fxq "$genome" $REFDIR/salinilactam_search/confirmed_genomes.slm.txt
	then
		cp $REFDIR/input_files/${cluster}* $REFDIR/salinilactam_search/tempffn/
		genbank_to_fasta.py -i ${cluster}.gbk -m genbank -o ${cluster}.ffn -s nt > /dev/null 2>&1
	else
		:
	fi

done < $REFDIR/salinilactam_search/slm_update.txt

# combine files from same genome
ls S*.ffn | cut -f1 -d'.' | sort -u > genomes.txt

while read genome 
do
	genomename=$(echo $genome | sed 's/_//g')
	# add in genome name before it gets lost
	cat ${genome}.cluster*.ffn | awk -v var="$genomename" '/^>/{sub(">", ">"var"_"++i " ")}1' > ${genome}.fix.ffn
	rm ${genome}.cluster*.ffn
done < genomes.txt

# now start the clustering step
cd $REFDIR/salinilactam_search/confirmed_clusters2
bgccount=$(ls S*.fa | cut -f1 -d'.' | wc -l | tr -d "[:blank:]")

cd $REFDIR/salinilactam_search/tempffn/
cat *.fix.ffn > total.slm.fasta
rm -f *.fix.ffn
rm -f *.gbk

cd-hit-est -i total.slm.fasta -o slm -c 0.8 -g 1 -sc 1 > /dev/null 2>&1

clstr2txt.pl slm.clstr > slm.cdhit.txt

cat slm.cdhit.txt | awk -v "seq=$bgccount" '{if ($3 == seq) print $0;}' > slm.cdhit.allgenomes.txt

cat slm.cdhit.allgenomes.txt | cut -f2 | sort -u > slm.cdhit.allgenomes.clusters.txt

echo "aligning the genes found in all genomes in slm BGC"

while read line
do

	cat slm.cdhit.allgenomes.txt | awk -v "seq=$line" '{if ($2 == seq) print $0;}' | cut -f1 > slm.cdhit.allgenomes.cluster${line}.txt
	geneIDsearch=$(cat slm.cdhit.allgenomes.txt | awk -v "seq=$line" '{if ($2 == seq) print $0;}' | awk '{if ($5 ==1) print $1;}' )
	geneID=$(grep -w "$geneIDsearch" total.slm.fasta | cut -f2 -d'|' | awk '{$1=$1};1')
	genedescript=$(grep -w "$geneIDsearch" total.slm.fasta | cut -f3 -d'|' | awk '{$1=$1};1')

	# check to make sure each gene is only represented by one genome
	function countuniquelines () { sort "$1" | cut -f1 -d'_' | uniq -c | sort -bgr; }
	check=$(countuniquelines slm.cdhit.allgenomes.cluster${line}.txt | head -n1 | tr -s ' ' | cut -f2 -d' ' )

	if [ $check -ne "1" ]
	then

		echo -e "salinilactam\tcluster${line}\t${geneIDsearch}\t${geneID}" >> $REFDIR/coreBGCgenes/repetitivegenes.txt

	else
		search-fasta.py -i total.slm.fasta -m slm.cdhit.allgenomes.cluster${line}.txt -o total.slm.cluster${line}.temp.fasta

		sed '/^>/ s/_.*//' total.slm.cluster${line}.temp.fasta >> total.slm.cluster${line}.fasta
		clustalo -i total.slm.cluster${line}.fasta -o total.salinilactam.cluster${line}.aln --force

		count-bp-in-seqs.sh total.salinilactam.cluster${line}.aln > /dev/null 2>&1
		bpcount=$(grep -v ">" total.salinilactam.cluster${line}-bp-count.txt | grep -v "sequences" | sort | uniq)

				#	# check to make sure open reading frames are good to go!
		if (( $bpcount % 3 == 0 ))           # no need for brackets
			then
				mv total.salinilactam.cluster${line}.aln $REFDIR/coreBGCgenes/goodalignments
			else
				mv total.salinilactam.cluster${line}.aln $REFDIR/coreBGCgenes/badalignments
		fi

		echo -e "salinilactam\tcluster${line}\t${bpcount}\t${geneID}\t${genedescript}" >> $REFDIR/coreBGCgenes/goodcoreBGCgenes.txt
		rm -f total.slm.cluster${line}.fasta
		rm -f total.slm.cluster${line}.temp.fasta
	fi

	rm -f slm.cdhit.allgenomes.cluster${line}.txt
	rm -f total.salinilactam.cluster${line}-bp-count.txt

done < slm.cdhit.allgenomes.clusters.txt

rm -f slm.cdhit.allgenomes.clusters.txt
rm -f slm.cdhit.allgenomes.txt






