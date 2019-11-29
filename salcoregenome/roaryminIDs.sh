#!/bin/bash

BASE=/Volumes/Bennett_BACKUP/Research/jensen-metabolites/salinispora-genomics

GENOMEDIR=$BASE/prodigal_annotation
WORKDIR=$BASE/roary

cd $GENOMEDIR

# run PROKKA
for f in *.fna
do

	genomename=${f%.fna}


	# rename contigs with BBMAP for PROKKA accomodations 
	rename.sh in=$f out=$genomename.temp.fa prefix=contig > /dev/null 2>&1

	# translate using prokka, skip empty files
	if [[ ! -s ${genomename}.temp.fa  ]]
		then
			echo ${genomename} >> $GENOMEDIR/somethingwrong.txt
			continue
		else
			echo "Annotating ${genomename}..."

			prokka --kingdom Bacteria --outdir ${genomename} ${genomename}.temp.fa --quiet --genus Salinispora --prefix $genomename > /dev/null 2>&1

	fi

	rm $genomename.temp.fa


done


###	run ROARY on test genomes to get an idea for the right parameters to use

while read line
do

	cp $GENOMEDIR/$line/*.gff $WORKDIR/test_subsample

done < $WORKDIR/genome_subset.txt

find . -type f -name '*.gff' -exec cp {} $WORKDIR  \;

cd $WORKDIR

for minID in {85}
do

	rm -rf roary${minID}

	roary -f roary${minID} -p 8 -i ${minID} *.gff


	cd roary${minID}

	create_pan_genome_plots.R
	python $BASE/roary_plots.py accessory_binary_genes.fa.newick gene_presence_absence.csv 


done

rm *.gff


