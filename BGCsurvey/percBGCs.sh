#!/bin/bash

BASEDIR=/Volumes/Bennett_BACKUP/Research/jensen-metabolites/micromonospora_BGCsurvey/patric_genomes
REFDIR=$BASEDIR/antismash
GENOMEDIR=$BASEDIR/genomes

cd $REFDIR


source activate antismash

# set thread limit for parallel processes
# set the max background jobs to be run at a given time to not overkill processor
# you can set it to more than 

function max_bg_procs {
	if [[ $# -eq 0 ]] ; then
			echo "Usage: max_bg_procs NUM_PROCS.  Will wait until the number of background (&)"
			echo "  bash processes (as determined by 'jobs -pr') falls below NUM_PROCS"
			return
	fi
	local max_number=$((0 + ${1:-0}))
	while true; do
			local current_number=$(jobs -pr | wc -l)
			if [[ $current_number -lt $max_number ]]; then
					break
			fi
			sleep 1
	done
}

for f in *.gbk
do

	max_bg_procs 4
	genomename=${f%.gbk}

	# no need to do pairwise if reciprocal ANI values were aleady calculated
	if [ -d "$genomename" ]; then
		continue

	else
		antismash $f &>/dev/null; fi & done 

wait 

source deactivate antismash

MOORE=/Volumes/Bennett_BACKUP/Research/jensen-metabolites/micromonospora_BGCsurvey/moorea_subset/antismash
cd $MOORE

ls *.gbk | sed 's/.gbk//g' > genomes.txt

# echo -e "genomeID\tgenus\tgenomelength\ttotalBGClength\tpercBGCs\tnumBGCs" > $BASEDIR/percBGCs.txt
echo -e "genomeID\tgenus\tgenomelength\ttotalBGClength\tpercBGCs\tnumBGCs" > $MOORE/percBGCs.txt

while read genomename
do
	echo $genomename
	genus=$(echo $genomename | cut -f1 -d '_')

	cd ../${genomename}
	# cd $GENOMEDIR/${genomename}

	n50_calc.py $genomename.fna
	genomelength=$(cat ${genomename}_n50.txt | cut -f3 | grep -v "Genome")

	rm -f ${genomename}_n50.txt

	# cd $REFDIR/$genomename
	cd ../antismash/$genomename


	count=`ls -1 *.cluster*.gbk 2>/dev/null | wc -l`

	if [ $count != 0 ]
	then 

		for f in *.cluster*.gbk
		do

			gb2fasta.py $f fasta

		done

		cat *.cluster*.fasta > $genomename.total.fasta
		n50_calc.py $genomename.total.fasta
		

		numBGC=$(cat $genomename.total.fasta | grep '>' | wc -l | tr -d "[:blank:]" )
		BGClength=$(cat $genomename.total_n50.txt | cut -f3 | grep -v "Genome")

		calc(){ awk "BEGIN { print "$*" }"; }
		percBGC=$(calc $BGClength/$genomelength)

		rm -f *.cluster*.fasta
		rm -f $genomename.total.fasta
		rm -f $genomename.total_n50.txt
		
	else
		numBGC=$count
		BGClength=$count
		percBGC=$count
	fi 
	

	echo -e "${genomename}\t${genus}\t${genomelength}\t${BGClength}\t${percBGC}\t${numBGC}" >> $MOORE/percBGCs.txt

	# cd $REFDIR
	cd $MOORE

done < genomes.txt
