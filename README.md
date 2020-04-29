Relevant files for the analysis of the distribution of biosynthetic gene clusters (BGCs) in the marine sediment  bacterium, *Salinispora*

# salBGCevol
Evolution of the diversity and distribution of BGCs in the marine sediment bacterium Salinispora. Currently, there are 118 publicly available Salinispora genomes collected from around the world in predominantly marine sediments in tropical and sub-tropical areas. This project examines how this genus differentiated into the nine described species.

<p align="center">
  <img width="860" height="400" src="images/updated-species-map-01.png">
</p>

# Order of Analysis
1. Core and Flexible Genome

⋅⋅⋅Genomes were processed to find core and flexible genomes across the 118 strains. Analysis of the flexible genome revealed species-specific genes were clearly linked to biosynthetic potential to produce secondary metabolites.

2. Biosynthetic Gene Clusters (BGCs) - antiSMASH

⋅⋅⋅Strains were ran through antiSMASH to identify BGCs across all strains. The average number of BGCs per Salinispora strain was 25.8 BGC fragments.

3. BGC Survey

⋅⋅⋅Does Salinispora have elevated abundances of BGCs? To answer this, we compared the number of BGCs identified across a wide range of well-characterized bacterial genera (N=589 genomes from 41 genera). We standardized by genome size by dividing the total number of bp in identified BGCs (via antiSMASH v5) by the genome length, providing a metric we called "Percent of Genome Dedicated to Seconary Metabolism"

4. Evidence for HGT in BGCs

⋅⋅⋅We first compared the %GC and codon usage biases for all BGCs in the above analysis to Salinispora. We found little evidence to support recent inter-genus transfer of BGCs. Further, when we compared within Salinispora, we found rare instances of inter-species transfer (as defined by %ANI).

5. Phylogenetic conservatism of BGCs

⋅⋅⋅Due to the large biosynthetic potential within Salinispora, we needed to dereplicate the identified BGCs into gene cluster families (GCFs) using BiG-SCAPE. We identified 305 GCFs across the strains and the distribution of these GCFs exhibited a strong species signature.

6. Evolutionary processes contributing to BGC diversity

⋅⋅⋅To evaluate evolutionary processes, we focused on 9 BGCs where the biosynthesis and product were previously identified. We evaluated the distribution of these 9 BGCs and looked for evidence of recombination, selective sweeps, directional selection, and BGC pathway evolution.

7. Linking biosynthetic potential to metabolite production

⋅⋅⋅We grew 30 Salinispora strains in triplicate and extracted metabolites produced over a week time. We compared the full metabolome across strains with data processing in MZMine and GNPS. Further, we applied targeted metabolomics to extract the 9 BGCs analyzed above. We hypothesized that despite the presence of the BGC across species, subtle differences in the BGCs would contribute to differential production of analogs of the molecule.

# Software and Databases Used
[NCBI Salinispora Genomes](https://www.ncbi.nlm.nih.gov/genome/?term=salinispora)

## Core and Flexible Genome
Ortholog Prediction: [ROARY](https://sanger-pathogens.github.io/Roary/)

Phylogenetic Analysis: [RAxML](https://cme.h-its.org/exelixis/software.html)

## BGC Analysis
BGC Identification: [antiSMASH](https://antismash.secondarymetabolites.org)

> Originally ran with antiSMASH v4 and BiG-SCAPE v1 but since been updated to antiSMASH v5 and BiG-SCAPE v1.1

BGC Clustering: [BiG-SCAPE](https://git.wageningenur.nl/medema-group/BiG-SCAPE)

BGC Alignment: [Mugsy](http://mugsy.sourceforge.net/)

## Recombination Analyses
[ClonalFrameML](https://github.com/xavierdidelot/ClonalFrameML)
