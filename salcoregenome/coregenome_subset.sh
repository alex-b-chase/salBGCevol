#!/bin/bash


usage(){
	echo "$(basename "$0") [-h] 

This program to take output from ROARY and subset core identified core genomes from outdirXX


In other words, you cannot have .fa files that represent both nucleotide genomes and translated amino acid genomes

correct use:
	coregenome_subset.sh -i minID 

where:
	--help or -h  			show this help text
	--identity or -i 		min BLASTp ID that ROARY generated (output files must be named 'outdir{minID}')
"
}


# get user defined options
if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit 1
fi


# get user input for other parameters
while test $# -gt 0; do
		case "$1" in
				-i)
						shift
						if test $# -gt 0; then
								export IDENTITY=$1
						else
								echo "no identity specified"
								exit 1
						fi
						shift
						;;
				--identity*)
						export IDENTITY=`echo $1 | sed -e 's/^[^=]*=//g'`
						shift
						;;
				*)
						break
						;;
		esac
done




### going to subset core genome genes from all Salinispora

BASEDIR=/Volumes/Bennett_BACKUP/Research/jensen-metabolites/salinispora-genomics

### use the XX% BLASTp minID based on the results in: 
### tested on a subset of genomes and picked the BLAST pID

# get number of genomes for all
numgen=$(ls $BASEDIR/roary/*.gff | wc -l)

minID=$IDENTITY 

rm -f $BASEDIR/coregenes${minID}.txt

################################################################################################
# go through the core genes output from ROARY
# create reference files for each core gene across all genomes
# output the results into a gene-specific reference file
################################################################################################

csv2txt.py < $BASEDIR/roary/roary${minID}/gene_presence_absence.csv > $BASEDIR/coregenes${minID}.temp.txt

head -n1 $BASEDIR/coregenes${minID}.temp.txt | cut -f1,3,4,15- | tr -d \'\" > $BASEDIR/coregenes${minID}.txt

awk -F"\t" -v var="$numgen" '(NR>1) && ($4 == var ) ' $BASEDIR/coregenes${minID}.temp.txt | \
cut -f1,3,4,15- | tr -d \'\" | sed 's/ /_/g' >> $BASEDIR/coregenes${minID}.txt

cut -f1,2 $BASEDIR/coregenes${minID}.txt > $BASEDIR/coregeneIDs${minID}.txt

rm -f $BASEDIR/coregenes${minID}.temp.txt


# transpose the file
cut -f4- $BASEDIR/coregenes${minID}.txt | awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' | tr -d '\r' | head -n${numgen} > $BASEDIR/coregenes${minID}.ids.txt

rm -rf $BASEDIR/tmp
mkdir -p $BASEDIR/tmp

while read genome genes
do

	echo $genes | tr ' ' '\n' > $BASEDIR/tmp/$genome.txt

done < $BASEDIR/coregenes${minID}.ids.txt

# rm -f $BASEDIR/coregenes${minID}.ids.txt

cd $BASEDIR/tmp/

for f in *.txt
do

	genome=${f%.txt}

	echo -ne "starting with ${genome}"

	while read geneID
	do

		gene=$(grep -w ${geneID} $BASEDIR/coregenes${minID}.txt | cut -f1 | sed 's;/;_;g')

		if [[ $geneID =~ [A-Z] ]]
			then

			locusID=${geneID}

			else
				locusID=$(grep "ID=${geneID};" $BASEDIR/roary/roary${minID}/fixed_input_files/${genome}.gff | \
				grep -o 'locus_tag=.*' | cut -f2 -d'=' | cut -f1 -d';')

		fi

		echo -e "${genome}\t${locusID}" >> ${gene}.ids.txt

	done < $f
	

	echo -e "\t - done!"

done


echo ""
echo "Starting the organization by each core gene"
echo ""


################################################################################################
# go through the reference ids.txt files to find the genome file and sequence corresponding
# to the AA sequence of the core genes
# output the results into a core gene specific AA file with all seqs in there
################################################################################################

# rm -rf $BASEDIR/coregenes
mkdir -p $BASEDIR/coregenes
mkdir -p $BASEDIR/coregenes/goodcore
mkdir -p $BASEDIR/coregenes/geneIDs
mkdir -p $BASEDIR/coregenes/nucleotide

# mv $BASEDIR/tmp/*.ids.txt $BASEDIR/coregenes
rm -rf $BASEDIR/tmp/

cd $BASEDIR/coregenes

# remove old files before beginning
rm -f *.total.aln
rm -f *.total.fa 
rm -f *.fasta 
rm -f *.faa 
rm -f *.fas
rm -f *.ffn 
rm -f *.aln

# loop through the lines with the hit core genes and subset by the specific core gene
for f in *.ids.txt
do
	filename=${f%.ids.txt}

	echo -ne "Processing core gene ${filename}"

	# check to make sure these are not paralogs
	check=$(cut -f1 $f | sort | uniq -c | grep -v '^ *1 ')

	if [ -z "$check" ]
		then
			:
		else 
			echo "$filename is a paralog. Deleting..."
			echo $filename >> $BASEDIR/paralogs.txt 

			rm $f 
			break
	fi

	echo -ne "...subsetting"

	while read genome sequence 
	do
		genfile=${genome}.faa
		genfile2=${genome}.ffn
		sequence2=$(echo $sequence | cut -f1 -d' ')

		# some genomes may not exist for some reason or were not included 
		# so if genome was deleted, do nothing (:)
		if [ -f $BASEDIR/prodigal_annotation/$genome/$genfile2 ]
			then 
				# awk -v "seq=$sequence2" -v RS='>' '$1 == seq {print RS $0}' $BASEDIR/prodigal_annotation/$genome/$genfile | \
				# sed "s/[>].*$/>$genome/" >> $BASEDIR/coregenes/$filename.ids.faa

				awk -v "seq=$sequence2" -v RS='>' '$1 == seq {print RS $0}' $BASEDIR/prodigal_annotation/$genome/$genfile2 | \
				sed "s/[>].*$/>$genome/" >> $BASEDIR/coregenes/$filename.ids.ffn

				# echo ${sequence2} > $BASEDIR/coregenes/$filename.$genome.temp.txt

				# filterbyname.sh in=$BASEDIR/prodigal_annotation/$genome/$genfile2 \
				# out=$BASEDIR/coregenes/$filename.$genome.temp.ffn \
				# include=t names=$BASEDIR/coregenes/$filename.$genome.temp.txt \
				# substring=name ow=t > /dev/null 2>&1

				# sed "s/[>].*$/>$genome/" $BASEDIR/coregenes/$filename.$genome.temp.ffn >> $BASEDIR/coregenes/$filename.ids.ffn

			else
				:
		fi

		rm -f $BASEDIR/coregenes/$filename.$genome.temp.ffn
		rm -f $BASEDIR/coregenes/$filename.$genome.temp.txt

	done < $f

	# cat $BASEDIR/coregenes/$filename.ids.faa | tr -d "[ -%,;\(\):=\.\\\[]\"\']" | sed "s/\*//g" > $BASEDIR/coregenes/$filename.total.fa 
	cat $BASEDIR/coregenes/$filename.ids.ffn | tr -d "[ -%,;\(\):=\.\\\[]\"\']" | sed "s/\*//g" > $BASEDIR/coregenes/$filename.total.fas 

	echo -ne "...aligning"

	# align each gene
	# clustalo -i $BASEDIR/coregenes/$filename.total.fa -o $BASEDIR/coregenes/$filename.total.faa.aln --force
	clustalo -i $BASEDIR/coregenes/$filename.total.fas -o $BASEDIR/coregenes/$filename.total.ffn.aln --force

	echo -ne "...checking\t"

	check2=$(grep '>' $BASEDIR/coregenes/$filename.total.ffn.aln | wc -l)
	if [[ $check2 == $numgen ]]
		then
			mv $BASEDIR/coregenes/$filename.total.ffn.aln $BASEDIR/coregenes/goodcore
			mv $f $BASEDIR/coregenes/geneIDs
			mv $BASEDIR/coregenes/$filename.total.fas $BASEDIR/coregenes/nucleotide
			echo $filename >> goodcoregenes.txt
		else 
			no_dups.sh $BASEDIR/coregenes/$filename.total.ffn.aln
			echo "check on ${filename}.total.aln"
	fi

	rm -f $BASEDIR/coregenes/$filename.ids.faa
	rm -f $BASEDIR/coregenes/$filename.ids.ffn
	echo -e "- done!"

done



