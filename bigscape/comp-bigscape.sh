#!/bin/bash

REFDIR=/Volumes/Bennett_BACKUP/Research/jensen-metabolites/salinispora-genomics/antismash5
OUTDIR=$REFDIR/bigscape/input_files

mkdir -p $OUTDIR


while read genome
do

	echo "processing $genome ..."

	cd $REFDIR/$genome 

	for f in contig*.gbk 
	do
		cp $f $OUTDIR/$genome.$f
	done

done < $REFDIR/genomes.txt

cd $REFDIR/bigscape
# rm -rf $REFDIR/bigscape/bigscape_out

source activate bigscape

BIGDIR=/anaconda2/envs/bigscape/BiG-SCAPE
PFAMDIR=/Users/alexchase/software/PfamScan

python $BIGDIR/bigscape.py \
-i input_files \
-o bigscape_out \
--pfam_dir $PFAMDIR \
-c 6 --mibig --cutoffs 0.3 \
--clan_cutoff 0.3 0.7 \
--mix --no_classify

# can try with different cutoff thresholds to get BGC families in order

# for cutoffvaluemin in {0.2, 0.25, 0.3, 0.35, 0.4}
# do

# 	python $BIGDIR/bigscape.py \
# 	-i input_files \
# 	-o bigscape_out \
# 	--pfam_dir $PFAMDIR \
# 	-c 6 --mibig --cutoffs 0.15 \
# 	--clan_cutoff $cutoffvaluemin 0.7 \
# 	--include_singletons --verbose

# done

source deactivate bigscape
