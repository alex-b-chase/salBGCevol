## Identification of biosynthetic gene clusters

Annotated genomes from PROKKA were passed through antiSMASH v5 for the identification of BGCs. 

The resulting output was fed to BiG-SCAPE to dereplicate BGCs into gene cluster families (GCFs) and match with known BGCs in the curated MIBIG database.

BiG-SCAPE computes a similarity score between each pairwise BGC. This data can be used to visualize the network:

<p align="center">
  <img width="460" height="500" src="../images/bgcnetwork-01.png">
</p>
