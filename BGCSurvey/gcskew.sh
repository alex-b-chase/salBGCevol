#!/bin/bash

EXECDIR=/Users/alexchase/software/codonW 

BASE=/Volumes/Bennett_BACKUP/Research/jensen-metabolites/salinispora-BGCevol
BASEDIR=$BASE/genomic-data
REFDIR=$BASE/antismash
GENOMEDIR=$BASEDIR/prokka_annotation

OUTDIR=$BASE/BGCsurvey
metadata=$BASE/genomeID.txt

cd $REFDIR

# echo -e "genome\tgenus\tclade\ttype\tregion\tGC%" > $OUTDIR/gcskew.txt
# echo -e "genome\tgenus\tclade\ttype\tregion\tCBImean\tCAImean\tFopmean" > $OUTDIR/codonbias.txt

# # do Salinispora first
# while read genome
# do

# 	cd $GENOMEDIR/$genome
# 	n50_calc.py ${genome}.fna
# 	gcgen=$(cut -f6 ${genome}_n50.txt | grep -v "GC")
# 	clade=$(grep -w "${genome}" $metadata | cut -f5)

# 	# get average genome CBI value
# 	$EXECDIR/codonw ${genome}.ffn -all_indices -nomenu -silent > /dev/null 2>&1
# 	cbigenome=$(grep -v "title" ${genome}.out | awk '{ sum += $7; n++ } END { if (n > 0) print sum / n; }')
# 	caigenome=$(grep -v "title" ${genome}.out | awk '{ sum += $6; n++ } END { if (n > 0) print sum / n; }')
# 	fopgenome=$(grep -v "title" ${genome}.out | awk '{ sum += $8; n++ } END { if (n > 0) print sum / n; }')

# 	echo -e "${genome}\tSalinispora\t${clade}\twholegenome\twholegenome\t${gcgen}" >> $OUTDIR/gcskew.txt
# 	echo -e "${genome}\tSalinispora\t${clade}\twholegenome\twholegenome\t${cbigenome}\t${caigenome}\t${fopgenome}" >> $OUTDIR/codonbias.txt

# 	rm -f ${genome}.out
# 	rm -f ${genome}.blk
# 	rm -f ${genome}_n50.txt

# 	cd $REFDIR/${genome}
# 	for bgccluster in *.region*.gbk
# 	do
# 		bgccluster2=${bgccluster%.gbk}
# 		genbank_to_fasta.py -i $bgccluster -m genbank -o ${bgccluster2}.fasta -s whole > /dev/null 2>&1
# 		n50_calc.py ${bgccluster2}.fasta

# 		gcBGC=$(cut -f6 ${bgccluster2}_n50.txt | grep -v "GC")
# 		echo -e "${genome}\tSalinispora\t${clade}\tBGC\t${bgccluster2}\t${gcBGC}" >> $OUTDIR/gcskew.txt

# 		# get BGC codon data
# 		genbank_to_fasta.py -i $bgccluster -m genbank -o ${bgccluster2}.ffn -s nt > /dev/null 2>&1
# 		$EXECDIR/codonw ${bgccluster2}.ffn -all_indices -nomenu -silent > /dev/null 2>&1

# 		cbiBGC=$(grep -v "title" ${bgccluster2}.out | awk '{ sum += $7; n++ } END { if (n > 0) print sum / n; }')
# 		caiBGC=$(grep -v "title" ${bgccluster2}.out | awk '{ sum += $6; n++ } END { if (n > 0) print sum / n; }')
# 		fopBGC=$(grep -v "title" ${bgccluster2}.out | awk '{ sum += $8; n++ } END { if (n > 0) print sum / n; }')

# 		echo -e "${genome}\tSalinispora\t${clade}\tBGC\t${bgccluster2}\t${cbiBGC}\t${caiBGC}\t${fopBGC}" >> $OUTDIR/codonbias.txt

# 		rm -f ${bgccluster2}.ffn
# 		rm -f ${bgccluster2}.out
# 		rm -f ${bgccluster2}.blk
# 		rm -f ${bgccluster2}.fasta
# 		rm -f ${bgccluster2}_n50.txt
# 	done

# done < $REFDIR/genomes.txt

# grep -v "^SaCN" $OUTDIR/patric_genomes/antismash/genomes.txt | grep -v "^SaDSM" | grep -v "^SpCN" | \
# grep -v "^StCN" > $OUTDIR/patric_genomes/antismash/genomes.temp.txt

# do the BGC survey genomes next
while read genome
do

	cd $OUTDIR/patric_genomes/prokka_annotation/$genome
	n50_calc.py ${genome}.fna
	gcgen=$(cut -f6 ${genome}_n50.txt | grep -v "GC")
	genus=$(echo $genome | cut -f1 -d'_')

	# get average genome CBI value
	$EXECDIR/codonw ${genome}.ffn -all_indices -nomenu -silent > /dev/null 2>&1
	cbigenome=$(grep -v "title" ${genome}.out | awk '{ sum += $7; n++ } END { if (n > 0) print sum / n; }')
	caigenome=$(grep -v "title" ${genome}.out | awk '{ sum += $6; n++ } END { if (n > 0) print sum / n; }')
	fopgenome=$(grep -v "title" ${genome}.out | awk '{ sum += $8; n++ } END { if (n > 0) print sum / n; }')

	echo -e "${genome}\t${genus}\tNA\twholegenome\twholegenome\t${gcgen}" >> $OUTDIR/gcskew.txt
	echo -e "${genome}\t${genus}\tNA\twholegenome\twholegenome\t${cbigenome}\t${caigenome}\t${fopgenome}" >> $OUTDIR/codonbias.txt

	rm -f ${genome}.out
	rm -f ${genome}.blk
	rm -f ${genome}_n50.txt

	cd $OUTDIR/patric_genomes/antismash/$genome

	count=`ls -1 *.region*.gbk 2>/dev/null | wc -l`

	if [ $count != 0 ]
	then 
		for bgccluster in *.region*.gbk
		do
			bgccluster2=${bgccluster%.gbk}
			genbank_to_fasta.py -i $bgccluster -m genbank -o ${bgccluster2}.fasta -s whole > /dev/null 2>&1
			n50_calc.py ${bgccluster2}.fasta

			gcBGC=$(cut -f6 ${bgccluster2}_n50.txt | grep -v "GC")
			echo -e "${genome}\t${genus}\tNA\tBGC\t${bgccluster2}\t${gcBGC}" >> $OUTDIR/gcskew.txt

			# get BGC codon data
			genbank_to_fasta.py -i $bgccluster -m genbank -o ${bgccluster2}.ffn -s nt > /dev/null 2>&1
			$EXECDIR/codonw ${bgccluster2}.ffn -all_indices -nomenu -silent > /dev/null 2>&1

			cbiBGC=$(grep -v "title" ${bgccluster2}.out | awk '{ sum += $7; n++ } END { if (n > 0) print sum / n; }')
			caiBGC=$(grep -v "title" ${bgccluster2}.out | awk '{ sum += $6; n++ } END { if (n > 0) print sum / n; }')
			fopBGC=$(grep -v "title" ${bgccluster2}.out | awk '{ sum += $8; n++ } END { if (n > 0) print sum / n; }')

			echo -e "${genome}\t${genus}\tNA\tBGC\t${bgccluster2}\t${cbiBGC}\t${caiBGC}\t${fopBGC}" >> $OUTDIR/codonbias.txt

			rm -f ${bgccluster2}.ffn
			rm -f ${bgccluster2}.out
			rm -f ${bgccluster2}.blk
			rm -f ${bgccluster2}.fasta
			rm -f ${bgccluster2}_n50.txt
		done
	else
		:
	fi 

done < $OUTDIR/patric_genomes/antismash/genomes.temp.txt



