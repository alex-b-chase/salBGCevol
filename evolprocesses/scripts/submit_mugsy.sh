#!/bin/bash


### first need to take the confirmed clusters from the BiG-SCAPE output and convert to full length fasta files
### genbank_to_fasta.py -i $f -m genbank -o ${f%.gbk}.fasta -s whole > /dev/null 2>&1

REFDIR=/dfs5/bio/abchase/BGCs
OUTDIR=$REFDIR/mugsy_out

cd $REFDIR

rm -rf $OUTDIR
mkdir -p $OUTDIR

while read BGC BGCfullname
do

		echo "#!/bin/bash
#$ -N ${BGC}_MUGSY
#$ -m ae
#$ -q abio128
#$ -ckpt blcr
#$ -pe openmp 2

module load perl/5.26.1
module load gcc/6.4.0

source /dfs5/bio/abchase/mugsy_x86-64-v1r2.3/mugsyenv.sh

REFDIR=/dfs5/bio/abchase/BGCs
OUTDIR=/dfs5/bio/abchase/BGCs/mugsy_out/${BGC}

rm -rf \$OUTDIR
rm -f \$REFDIR/mugsy_out/${BGC}.aln
mkdir -p \$OUTDIR
cd \$REFDIR/${BGC}

mugsy --directory \${OUTDIR} --prefix ${BGC} *.fasta

cd \$OUTDIR
#######Extract an alignment, e.g., the first LCB (locally co-linear block)
maf2fasta.pl 1 < ${BGC}.maf > ${BGC}.fasta

#######now need to remove \"= sign\" from end of files 
grep -v \"=\" ${BGC}.fasta > \$REFDIR/mugsy_out/${BGC}.aln

	
	" > ${BGC}_MUGSY.sh

	# qsub ${BGC}_MUGSY.sh
	# sleep 5s

done < BGCnames.txt
