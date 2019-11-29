#!/bin/bash

BASE=/Volumes/Bennett_BACKUP/Research/jensen-metabolites/micromonospora_BGCsurvey
OUTDIR=$BASE/patric_genomes
REFDB=/Volumes/Bennett_BACKUP/Research/reference_db
ESSEN=$REFDB/essen_singlecopy

mkdir -p $OUTDIR
cd $OUTDIR

rm -rf $OUTDIR/prokka_annotation
# rm -rf $OUTDIR/markergenes
# rm -rf $OUTDIR/genomes
rm -rf $OUTDIR/ftp.patricbrc.org
# rm -rf $OUTDIR/antismash
# mkdir -p $OUTDIR/genomes
mkdir -p $OUTDIR/prokka_annotation
# mkdir -p $OUTDIR/markergenes
# mkdir -p $OUTDIR/antismash

cut -f1 $BASE/PATRIC_genome.metadata.txt | grep -v "Genome" > $OUTDIR/download.txt

echo "downloading genomes from PATRIC"
# download the files from PATRIC
while read line
do
	filename=$line.fna

	cd $OUTDIR/genomes

	# if we already downloaded the file, skip it
	if [ ! -f "$filename" ]
	then
		cd $OUTDIR
		wget -cqNr -t 45 -A "*.fna" "ftp://ftp.patricbrc.org/genomes/${line}/" -P .
		# echo "done with $line"

	else
		# echo "done with $line"
		continue
	fi

done < download.txt

# rm -f download.txt
cd $OUTDIR/ftp.patricbrc.org

# move downloaded genomes to new directory
find . -type f -name '*.fna' -exec mv {} $OUTDIR/genomes \;

cd $BASE
rm -rf $OUTDIR/ftp.patricbrc.org

while read line
do
	cp -r  /Volumes/Bennett_BACKUP/Research/jensen-metabolites/salinispora-genomics/prodigal_annotation/${line} $OUTDIR/genomes/
	cp $OUTDIR/genomes/${line}/${line}.fna $OUTDIR/genomes/
done < salinispora.genomes.txt


# #################################################################################
# #################################################################################

cd $OUTDIR/genomes
# annotate and screen for marker genes
refmgs=$ESSEN/Essen_SingleCopy.txt
cat $refmgs | awk '{if ($6 == "1") print $0;}' | sed 's/ /_/g' > $OUTDIR/prokka_annotation/markers.sub.txt

mgs=$OUTDIR/prokka_annotation/markers.sub.txt
check=$(cat $mgs | wc -l)


# subset for naming purposes only
cut -f1-2 $BASE/PATRIC_genome.metadata.txt | sort -u > $OUTDIR/patric_list.txt
genomeMAP=$OUTDIR/patric_list.txt



# remove old files
rm -f *temp*
rm -f $OUTDIR/deleted_genomes.txt
rm -f $OUTDIR/excess_markers.txt
rm -f $OUTDIR/duplicate_markers.txt
rm -f $OUTDIR/final_passed_genomes.txt
rm -f $OUTDIR/somethingwrong.txt
rm -f $OUTDIR/genomes/*.markers.faa
rm -f $OUTDIR/genomes/*.faa
rm -f *ids.txt
rm -rf $OUTDIR/genomes/genome_markers
mkdir -p $OUTDIR/genomes/genome_markers

sleep 2s

genomecount=$(find . -maxdepth 1 -type f -name '*.fna'  | wc -l)
echo "DONE! There are ${genomecount} genomes ready to be annotated!"

# use PATRIC annotation and screen for marker genes!

for f in *.fna
do
	initial_start=$(date +%s)

	base=${f%.fna}
	patric=$(echo $f | cut -f1 -d'.')
	# check if genome is PATRIC or not - PATRIC genomes start with numbers
	if ! [[ $patric =~ ^[0-9]+$ ]]
	then
		output=$base
	else
		output=$(cat $genomeMAP | grep -w "$base" | cut -f2 | tr -d "[-/%,;\(\):=\.\\\[]\"\']" | sed 's/  */_/g')

		# if we can't find the $output variable in the mappingFile, output that for reference - should be nothing!
		if [ -z "$output" ] 
		then
			echo $base
			echo $output
			echo ''
		else
			:
		fi
	fi

	echo $output
	prokka \
	--outdir ${output} --prefix ${output} \
	--genus $(echo $output | cut -f1 -d'_') \
	--cpus 8 --fast --quiet $f

	cp ${output}/${output}.faa ./ 

	# search against the Banfield marker gene database
	hmmsearch --tblout $output.temp.txt -E 1e-05 --cpu 6 $ESSEN/Essen_SingleCopy.hmm $output.faa > /dev/null

	# sort HMMer output by best bit score for each gene
	grep -v "#" $output.temp.txt | tr -s ' ' | sort -k3,3 -k6,6nr | sort -u -k3,3 > $output.temp2.txt

	# now we can subset by the phylomarkers
	while read acc bla blah blah blahh blahhh protein protnum
	do
		grep $acc $output.temp2.txt | cut -f1 -d' ' | sort -u >> $output.temp3.txt

		filename=$(echo $output.faa)
		sequence=$(cat $output.temp2.txt | grep $acc | cut -f1 -d' ' | sort | uniq)
		prot=$(echo $protein | tr -d "[-/%,;\(\):=\.\\\[]\"\']")
		echo -e $filename'\t'$sequence >> ${prot}.ids.txt

	done < $mgs

	# check to make sure genomes only have one copy of each...
	protcheck=$(cat $output.temp3.txt | wc -l)

	if [ $protcheck -ne $check ]
	then
		echo -e $f'\t'$output'\t'$protcheck >> $OUTDIR/deleted_genomes.txt
		rm -rf *${output}*
		continue

	else
		echo "Genome ${output} passed!"
		cat $BASE/PATRIC_genome.metadata.txt | grep -w "$base" >> $OUTDIR/final_passed_genomes.txt
		cp ${output}/${output}.gbk $OUTDIR/antismash/
		

	fi		

	# filter out the marker genes from the genome
	search-fasta.py -i $output.faa -m $output.temp3.txt -o $output.markers.faa
	mv $output.markers.faa $OUTDIR/genomes/genome_markers/
	mv $output.faa $OUTDIR/prokka_annotation/


	total_end=$(date +%s)
	total_runtime=$(echo "$total_end - $initial_start" | bc -l)

	echo "Done annotating $output - total time $total_runtime seconds"

	rm *.temp*


done

find . -type f -maxdepth 1 -name '*.ids.txt' -exec cp {} $OUTDIR/prokka_annotation \;

###############################################################################################
# go through the reference ids.txt files to find the genome file and sequence corresponding
# to the AA sequence of the marker genes
# output the results into a marker gene specific AA file with all seqs in there
###############################################################################################

cd $OUTDIR/prokka_annotation

genomecount=$(find . -type f -name '*.faa'  | wc -l)
echo "DONE! There are ${genomecount} genomes that were downloaded and processed from PATRIC"

# remove all old files that may have previously been ran
rm -f *.ids.faa


# loop through the lines with the hit marker genes and subset by the marker genes
for f in *.ids.txt
do
	filename=${f%.ids.txt}

	echo "Processing marker gene $filename"

	while read gen sequence 
	do
		genome=$(echo $gen | cut -f1 -d'.')
		sequence2=$(echo $sequence | cut -f1 -d' ')
		#echo $sequence2

		# some genomes were deleted 
		# so if genome was deleted, do nothing (:)
		if [ -f ${genome}.faa ]
			then 
				awk -v "seq=$sequence2" -v RS='>' '$1 == seq {print RS $0}' ${genome}.faa | \
				sed "s/[>].*$/>$genome/" >> $filename.ids.faa
			else
				:
		fi
	done < $f
done

# # remove all old folders and files related to this
# # move the newly produced AA files to a new folder

cd $OUTDIR/prokka_annotation

cp *.ids.faa $OUTDIR/markergenes
cp *.ids.txt $OUTDIR/markergenes

# ################################################################################################
# # align all the individual protein files for tree construction
# ################################################################################################


cd $OUTDIR/markergenes

rm -f *.aln
rm -f *-no-dups.fasta


for f in *.ids.faa
do
	newname="${f/%.ids.faa/.total.fa}"
	alignment="${f/%.ids.faa/.total.aln}"

	# remove all weird characters (RAXML does not like them!)
	cat $f | tr -d "[ -%,;\(\):=\.\\\[]\"\']" | sed "s/\*//g" > "$newname"

	# align the sequences to create aligned fasta files
	clustalo -i $newname -o $alignment

	# remove duplicate IDs just in case there were any relic mistakes
	no_dups.sh $alignment
	# this will output our final file == "*-no-dups.fasta"

	# rm -f $f
done



# ################################################################################################
# # build a consensus reference tree by combining all 30 marker genes into one file
# # construct the tree
# ################################################################################################

# combine the fasta files together for tree construction keeping the proteins correctly ordered
# all individual alignment files must have the same header

catfasta2phyml.pl -f *.fasta  > concat.aligned.filtered.fa

rm -rf *.fasta

# build the tree using the aligned file and output 100 bootstrap replications
# kept crashing, running on the other lab computer - WORKED so run separately!

rm -f RAxML*
rm -rf $OUTDIR/markergenes/master_tree/
mkdir $OUTDIR/markergenes/master_tree/

mv concat.aligned.filtered.fa $OUTDIR/markergenes/master_tree/
cd $OUTDIR/markergenes/master_tree/

# FastTree concat.aligned.filtered.fa > concat.aligned.filtered.tre 

# try and run on HPC, no way will this work on desktop
raxml -s concat.aligned.filtered.fa -m PROTGAMMABLOSUM62 -n markergene -x 100 -# 100 -p 4321 -f a -T 6








