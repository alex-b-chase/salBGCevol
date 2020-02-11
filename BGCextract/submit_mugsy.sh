#!/bin/bash

REFDIR=/dfs3/bio/abchase/BGCs
OUTDIR=$REFDIR/mugsy_out

cd $REFDIR

rm -rf $OUTDIR
mkdir -p $OUTDIR

while read BGC BGCfullname
do

		echo "#!/bin/bash
#$ -N ${BGC}_MUGSY
#$ -m ae
#$ -q mic
#$ -pe openmp 2

module load perl/5.26.1
module load gcc/6.4.0

source /dfs3/bio/abchase/mugsy_x86-64-v1r2.3/mugsyenv.sh

REFDIR=/dfs3/bio/abchase/BGCs
OUTDIR=/dfs3/bio/abchase/BGCs/mugsy_out/${BGC}

mkdir -p \$OUTDIR
cd \$REFDIR/${BGC}

mugsy --directory \${OUTDIR} --prefix ${BGC} *.fasta

cd \$OUTDIR
#######Extract an alignment, e.g., the first LCB (locally co-linear block)
maf2fasta.pl 1 < ${BGC}.maf > ${BGC}.fasta

#######now need to remove \"= sign\" from end of files 
grep -v \"=\" ${BGC}.fasta > \$REFDIR/mugsy_out/${BGC}.aln

	
	" > ${BGC}_MUGSY.sh

	qsub ${BGC}_MUGSY.sh
	sleep 5s

done < BGCnames.txt
