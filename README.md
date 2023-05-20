# Taxonomy of the genus Serratia
Author:  
- Ekaterina Sukhinina  

Supervisors:  
- Kirill Antonets  
- Anton Shikov  


### Introduction  
In 2020 two new species of the genus Serratia were identified - Serratia nevei and S. bockelmannii. These species were categorized based on their biochemical and morphological characteristics. A phylogenetic tree was constructed using a limited set of genomes. It was found that many genomes annotated as S. marcescens likely belong to the two new species that were identified. These new species have not yet been officially recognized by taxonomists and are listed as not taxonomically validated in the NCBI database.  

#### Aim: 
To investigate whether the subdivision into species of the genus Serratia is consistent with the actual phylogeny

#### Objectives:
- Download the Serratia genomes from the NCBI Assembly database
- Construct a distance matrix between genomic assemblies using sourmash (doi: 10.21105/joss.00027)
- Build a table with the declared species and the closest species according to the distance matrix 
- Construct a phylogenetic tree based on the distance matrix.
- Obtain an alignment of the core genes (common to all assemblies). Based on this alignment, construct the phylogenetic tree.
- Determine at the genomic level the relevance of distinguishing new Serratis species. Determine the pangenome features of these species.  


### Workflow  
#### Uploading genomic assemblies
Get ids for all assemblies from ncbi_assembly by 'serratia' query:  
```
Entrez.email = "suhininaev@mail.ru"  
handle = Entrez.esearch(db="assembly", term='serratia', retmax='2600')  
record = Entrez.read(handle)  
ids = record['IdList']  
```
2449 ids (Date of access: 13.03.23)  

A python script (script.py) was written to download genomic assemblies from GenBank.

#### Reconstruction of a phylogenetic tree based on a distance matrix  
A distance matrix between genomic assemblies was constructed using sourmash (doi: 10.21105/joss.00027).  
Sourmash is a Python package that quickly compares potentially very large sets of DNA and protein sequences. This functionality can be used to, for example, cluster transcriptomes or genomes, to identify the taxonomy of new isolate or metagenome-assembled genomes, or to determine the taxonomic composition of a new metagenome sequence by comparing it against a database of reference genomes.  
Sourmash signatures contain one or multiple sub-sampled representations of DNA or protein sequences (FracMinHash sketches). Each FracMinHash sketch contains hashes (and optionally, hash abundances) that represent a subset of k-mers from the original sequences. The sourmash sketch command consistently subsamples k-mers across different sequences, so we can compare (e.g. intersect) sketches to understand sequence similarity. The command line function sourmash compare estimates similarity and containment metrics.
```
> /home/e_sukhinina/.local/bin/sourmash sketch dna -p scaled=1000,k=31 ../assemblies/*.fna.gz
> /home/e_sukhinina/.local/bin/sourmash compare -p 32 ./signatures/*.sig -o ./sourmash_results/serra_cmp --distance-matrix --ksize 31 --csv ./sourmash_results/dist_matrix.csv
```

The phylogenetic tree was reconstructed using a hierarchical clustering algorithm from the r stats package and ggtree package
```
dd <- as.dist(d)     # d - distance matrix
hclust_avg <- hclust(dd, method='average')
p1 <- ggtree(hclust_avg, layout="circular") 
plot(p1)
```

#### Phylogenetic tree annotation
The resulting phylogenetic tree was annotated according to the following parameters: attributed species, nearest species, taxonomy check status, type strain. We used ANI_report_prokaryotes.txt from the ncbi database to obtain this data:
https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/
ANI_report_prokaryotes.txt provides Average Nucleotide Identity (ANI) data that can be used to evaluate the taxonomic identity of genome assemblies of interest.  Also included is the ANI status which GenBank uses as a basis for decisions about taxonomic identity of public genome assemblies. Specific methods used can be found here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6978984/  
This file contains ANI data for all latest archaeal & bacterial genome assemblies in GenBank, together with the paired RefSeq genome assemblies when they exist.  
Accessing the data 5 April 2023.  
Reconstruction of an annotated tree. The annotation layers were superimposed manually in Adobe Photoshop, as there are too many genomes to build a clear heatmap.
```
p1 <- ggtree(hclust_avg, layout="circular") %<+% tax_table +
  geom_tippoint(aes(color=declared_species), size=6)+ 
  scale_color_manual(values = custom.col)
plot(p1)
```
The information on the type strains was taken from the NCBI taxonomy page: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Undef&id=613&lvl=3&keep=1&srchmode=1&unlock  
After clicking on the name of the species (species, not strain), we got a link to Genome in the table at the top right of the page. Clicking on it provided the Genbank ID of the assembly for the reference strain. It is important to note that the type strains are not available for all species. We have not considered unclassified Serratia, only the major species.

![alt text](https://github.com/shkot190819/serratia/blob/main/tree.jpg?raw=true)


