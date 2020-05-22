#!/bin/bash

BASE=/Volumes/Bennett_BACKUP/Research/jensen-metabolites/salinispora-BGCevol/BGCanalysis/evolprocesses

cd $BASE/mugsy_out

rm -f clonalframe-summary.txt
rm -f totalBGCxgenome.txt

echo -e "BGC\tnum_genomes\tBGC_alignment_length\tR-theta\tdelta\tnu\trmratio\tdeltanu" > clonalframe-summary.txt
echo -e "genomeID\tBGC" > totalBGCxgenome.txt

for f in *.aln
do

	BGC=${f%.aln}

	rm -rf ${BGC}_mugsy
	mkdir -p ${BGC}_mugsy

	echo -n "starting on ${BGC}..."

	genomenum=$(grep '>' $f | wc -l)
	grep '>' $f | sed 's/>//g' | cut -f1 -d' '> temp.txt
	cut -f1 -d'.' temp.txt > temp2.txt

	echo -n "processing ${genomenum} genomes..."

	paste temp.txt temp2.txt > newnames.txt

	fasta-rename.py $f newnames.txt ${BGC}.final.fa

	count-bp-in-seqs.sh ${BGC}.final.fa &> out

	bpcount=$(cat ${BGC}.final-bp-count.txt | grep -v ">" | grep -v "sequences" | sort -u)

	echo -n "BGC alignment is ${bpcount} bp..."
	echo -n "aligning with RAxML..."

	# if we already downloaded the file, skip it
	filename=RAxML_bipartitionsBranchLabels.${BGC}BGC
	if [ ! -f "$filename" ]
	then
		raxml -s ${BGC}.final.fa -m GTRGAMMA -n ${BGC}BGC -x 100 -# 100 -p 4321 -f a -T 7 
	else
		:
	fi

	cp RAxML_bipartitionsBranchLabels.${BGC}BGC $BASE/DTL_analysis

	# raxml -s ${BGC}.final.fa -m GTRGAMMA -n ${BGC}BGC -x 100 -# 100 -p 4321 -f a -T 7 &> out

	echo -n "computing recombination rates..."

	ClonalFrameML RAxML_bipartitionsBranchLabels.${BGC}BGC ${BGC}.final.fa ${BGC}BGC &> out
	Rscript /Users/alexchase/software/ClonalFrameML/src/cfml_results.R ${BGC}BGC &> out

	calc(){ awk "BEGIN { print "$*" }"; }
	recombination=$(cat ${BGC}BGC.em.txt | grep 'theta' | cut -f2 )
	delta=$(cat ${BGC}BGC.em.txt | grep 'delta' | cut -f2 )
	delta1=$(calc 1/$delta)
	nu=$(cat ${BGC}BGC.em.txt | grep 'nu' | cut -f2 )
	rmratio=$(calc $recombination*$delta1*$nu)
	deltanu=$(calc $delta1*$nu)
	echo -n "r/m ratio is $rmratio"

	# get a table of genome x BGCs
	cat ${BGC}.final.fa | grep '>' | sed 's/>//g' | \
	awk -v var="${BGC}" 'BEGIN{OFS=IFS="\t"} {print $0, var}' >> totalBGCxgenome.txt

	echo "...done with ${BGC}!"

	echo -e "${BGC}\t${genomenum}\t${bpcount}\t${recombination}\t${delta1}\t${nu}\t${rmratio}\t${deltanu}" >> clonalframe-summary.txt

	rm -f temp.txt
	rm -f temp2.txt
	rm -f newnames.txt
	rm -f ${BGC}.final-bp-count.txt
	rm -f RAxML_bestTree.${BGC}BGC
	rm -f RAxML_bipartitions.${BGC}BGC
	rm -f RAxML_bootstrap.${BGC}BGC
	rm -f RAxML_info.${BGC}BGC
	rm -f out
	rm -f ${BGC}.final.fa.reduced
	# rm -f $f

	mv ${BGC}BGC.em.txt ${BGC}_mugsy
	mv ${BGC}BGC.position_cross_reference.txt ${BGC}_mugsy
	mv ${BGC}BGC.labelled_tree.newick ${BGC}_mugsy
	mv ${BGC}BGC.ML_sequence.fasta ${BGC}_mugsy
	mv ${BGC}.final.fa ${BGC}_mugsy
	mv ${BGC}BGC.importation_status.txt ${BGC}_mugsy
	mv ${BGC}BGC.cfml.pdf ${BGC}_mugsy

done 


