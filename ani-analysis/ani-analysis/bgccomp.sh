#!/bin/bash
#$ -N BGC_comp
#$ -q mic
#$ -pe openmp 32


module load ruby/2.5.3 
module load blast/2.8.1

files=.fasta
FILES='*'$files

THREAD=32
ANIDIR=/data/users/abchase/enveomics/Scripts

REFDIR=/dfs5/bio/abchase/BGCs/fastafiles
TEMP=$REFDIR/tempfiles

rm -rf $TEMP
rm -f $REFDIR/total.txt; touch $REFDIR/total.txt
mkdir -p $TEMP
# move everything to new temporary directory
find . -maxdepth 1 -name "$FILES" -exec cp {} $TEMP \;

cd $TEMP

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

for f in $FILES
do
	# to get names of genomes without file extension
	gen1=$( basename $f $files )
	genome1=$(echo $f | cut -f1 -d'.')
	echo "Starting genome ${gen1} ..."

	for h in $FILES
	do
		genome2=$(echo $h | cut -f1 -d'.')
		max_bg_procs $THREAD

		gen2=$( basename $h $files )
		# no need to do pairwise if reciprocal ANI values were aleady calculated
		if grep -Fq "${gen2};${gen1}" ${REFDIR}/total.txt; then
			continue
		elif [ "$genome1" == "$genome2" ]; then
			continue

		else
			$ANIDIR/ani.rb -1 $f -2 $h -T ${gen1}'@'${gen2}.tab -q &>/dev/null; fi & done 

	wait
	# wait for parellization to finish for genome1 or limit from user defined threads

	# combine recently processed .tab files to master - limit number of files created overall...
	# this will combine all .tab files for pairwise comparisons for one genome at a time (for the loop)
	for i in *.tab
	do
		filename=$( basename $i .tab )
		gen1=$(echo $filename | cut -f1 -d'@')
		gen2=$(echo $filename | cut -f2 -d'@')
		ani=$(cat $i | cut -f1)
		echo ${gen1}';'${gen2}';'${ani} > ${gen1}.${gen2}.final.txt & done
	wait

	cat *.final.txt >> $REFDIR/total.txt
	

	# remove old files
	rm *.tab
	rm *.final.txt
	
done


tr ';' '\t' < $REFDIR/total.txt > $REFDIR/total2.txt
mv $REFDIR/total2.txt $REFDIR/total.txt



