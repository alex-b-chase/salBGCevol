#!/bin/bash

NOTDIR=/dfs5/bio/abchase/Notung-2.9.1.2
WORKDIR=/dfs5/bio/abchase/BGCs/DTL-analysis

LOSS=1.0
TRANS=3.0
DUP=2.0

cd $WORKDIR

rm -rf *.rooting*

# first need to find the optimum root for the phylogeny
for f in RAxML_bipartitionsBranchLabels.*BGC
do 
	BGC=$(echo $f | cut -f2 -d'.' | sed "s/BGC//g")

	java -jar $NOTDIR/Notung-2.9.1.2.jar -g $f \
	-s RAxML_boot.core.root.tre \
	--root --speciestag prefix \
	--costdup $DUP --costloss $LOSS \
	--maxtrees 10

	for (( i = 0; i <= 9; ++i ))
	do
	    echo "#!/bin/bash
#$ -N DTL-${BGC}.${i}
#$ -q bio
#$ -pe openmp 2

NOTDIR=$NOTDIR
WORKDIR=$WORKDIR

BGC=${BGC}

LOSS=1.0
TRANS=3.0
DUP=2.0

cd \$WORKDIR

java -jar \$NOTDIR/Notung-2.9.1.2.jar \\
--reconcile --speciestag prefix \\
--costdup \$DUP --costloss \$LOSS \\
--infertransfers true --costtrans \$TRANS \\
--parsable --treeoutput newick \\
-g ${f}.rooting.${i} \\
-s RAxML_boot.core.root.tre 

	    " > ${BGC}.${i}.DTL.sh
	    qsub ${BGC}.${i}.DTL.sh
	    
	done

done

