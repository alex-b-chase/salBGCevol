#!/bin/bash

BASEDIR=/Volumes/Bennett_BACKUP/Research/jensen-metabolites/salinispora-BGCevol/BGCsurvey/patric_genomes
REFDIR=$BASEDIR/antismash
GENOMEDIR=$BASEDIR/prokka_annotation

cd $REFDIR

export PATH="/anaconda2/bin:$PATH"
eval "$(conda shell.bash hook)"
conda activate antismash

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

find ./* -maxdepth 0 -type d | sed 's;./;;g' > genomes.txt

# for f in *.gbk
# do

# 	max_bg_procs 4
# 	genomename=${f%.gbk}

# 	# no need to do if already done
# 	if [ -d "$genomename" ]; then
# 		continue

# 	else
# 		antismash --genefinding-tool none $f &>/dev/null; fi & done 

# wait 

# conda deactivate 

echo -e "genomeID\tgenus\tgenomelength\ttotalBGClength\tpercBGCs\tnumBGCs" > $BASEDIR/percBGCs.txt

while read genomename
do
	echo $genomename
	genus=$(echo $genomename | cut -f1 -d '_')

	cd $GENOMEDIR/${genomename}

	n50_calc.py $genomename.fna
	genomelength=$(cat ${genomename}_n50.txt | cut -f3 | grep -v "Genome")

	rm -f ${genomename}_n50.txt

	cd $REFDIR/$genomename
	count=`ls -1 *.region*.gbk 2>/dev/null | wc -l`

	if [ $count != 0 ]
	then 

		for f in *.region*.gbk
		do

			gb2fasta.py $f fasta

		done

		cat *.region*.fasta > $genomename.total.fasta
		n50_calc.py $genomename.total.fasta
		

		numBGC=$(cat $genomename.total.fasta | grep '>' | wc -l | tr -d "[:blank:]" )
		BGClength=$(cat $genomename.total_n50.txt | cut -f3 | grep -v "Genome")

		calc(){ awk "BEGIN { print "$*" }"; }
		percBGC=$(calc $BGClength/$genomelength)

		rm -f *.region*.fasta
		rm -f $genomename.total.fasta
		rm -f $genomename.total_n50.txt
		
	else
		numBGC=$count
		BGClength=$count
		percBGC=$count
	fi 
	

	echo -e "${genomename}\t${genus}\t${genomelength}\t${BGClength}\t${percBGC}\t${numBGC}" >> $BASEDIR/percBGCs.txt

	cd $REFDIR

done < genomes.txt
