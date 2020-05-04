#!/bin/bash

BASEDIR=/Volumes/Bennett_BACKUP/Research/jensen-metabolites/salinispora-BGCevol
GENDIR=$BASEDIR/genomic-data/prokka_annotation
REFDIR=$BASEDIR/bigscape

cd $REFDIR

export PATH="/anaconda2/bin:$PATH"
eval "$(conda shell.bash hook)"

conda activate bigscape

BIGDIR=/anaconda2/envs/bigscape/BiG-SCAPE
PFAMDIR=/Users/alexchase/software/PfamScan

# while read BGC refMIBIG 
# do
# 	rm -rf $BASEDIR/BGCanalysis/individualBGCs/${BGC}/bigscape_out

# 	python $BIGDIR/bigscape.py -i input_files \
# 	-o $BASEDIR/BGCanalysis/individualBGCs/${BGC}/bigscape_out \
# 	--pfam_dir $PFAMDIR -c 6 --cutoffs 0.7 \
# 	--query_bgc $BASEDIR/BGCanalysis/individualBGCs/referenceBGCs/${refMIBIG}.gbk \
# 	--clan_cutoff 0.7 0.9 \
# 	--mix --no_classify --mode global

# done < bgc2refmibig.txt

BGC=salinilactam
refMIBIG=SpDSM45544_FJ.contig_0.region001

python $BIGDIR/bigscape.py -i input_files \
-o $BASEDIR/BGCanalysis/individualBGCs/${BGC}/bigscape_out \
--pfam_dir $PFAMDIR -c 6 --cutoffs 0.7 \
--query_bgc input_files/${refMIBIG}.gbk \
--clan_cutoff 0.7 0.9 \
--mix --no_classify --mode global
# --query_bgc $BASEDIR/BGCanalysis/individualBGCs/referenceBGCs/${refMIBIG}.gbk \


conda deactivate 

