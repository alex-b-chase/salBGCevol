## Genome Annotation

Salinispora genomes were first annotated with [PROKKA](https://github.com/tseemann/prokka) with a slight adjustment in the contig naming:

```bash
rename.sh in=$genome.fna out=$genome.temp.fa prefix=contig

prokka \
--outdir $BASEDIR/prokka-annotation/${genome} --prefix $genome \
--kingdom Bacteria --genus Salinispora \
$genome.temp.fa --quiet --force
```
## Core Genome Analysis

After annotation, we can call orthologous protein groups with [ROARY](https://github.com/sanger-pathogens/Roary). Feeding the output from PROKKA is very easy and compatible with ROARY. All we need is the resulting .gff files

However, we need to know the most crucial parameter, which is the minimum identity to cluster orthologs. We can use a subset of the genomes from all the species to figure this out.

```bash
for minID in {50,60,70,75,80,85,90,95}
do
	rm -rf roary${minID}
	roary -f roary${minID} -p 8 -i ${minID} *.gff

	cd roary${minID}
	create_pan_genome_plots.R
	# python $BASE/roary_plots.py accessory_binary_genes.fa.newick gene_presence_absence.csv 
done
```

Here is the output from this screening of minIDs:

<p align="center">
  <img width="460" height="300" src="https://github.com/alex-b-chase/salBGCevol/images/corepangenomeXblastID.png">
</p>

There was a nice inflection point at 85% so we can use all the genomes now and re-run ROARY with this parameter.

We can do some nice stats in R with the vegan package to look at how the presence/absence of flex genes structure Salinispora diversity:

<p align="center">
  <img width="460" height="500" src="https://github.com/alex-b-chase/salBGCevol/images/flexgenome.png">
</p>
