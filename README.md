# Taxonomy of the genus *Serratia*
Author:  
- Ekaterina Sukhinina  

Supervisors:  
- Kirill Antonets  
- Anton Shikov  


### Introduction  
In 2020 two new species of the genus *Serratia* were identified - *S. nevei* and *S. bockelmannii*. These species were categorized based on their biochemical and morphological characteristics. A phylogenetic tree was constructed using a limited set of genomes. It was found that many genomes annotated as *S. marcescens* likely belong to the two new species that were identified. These new species have not yet been officially recognized by taxonomists and are listed as not taxonomically validated in the NCBI database.  

#### Aim: 
To investigate whether the subdivision into species of the genus *Serratia* is consistent with the actual phylogeny

#### Objectives:
- Download the *Serratia* genomes from the NCBI Assembly database
- Construct a distance matrix between genomic assemblies using sourmash
- Build a table with the declared species and the closest bspecies according to the distance matrix 
- Construct a phylogenetic tree based on the distance matrix.
- Obtain an alignment of the core genes (common to all assemblies). Based on this alignment, construct the phylogenetic tree.
- Determine at the genomic level the relevance of distinguishing new Serratis species. Determine the pangenome features of these species.  


### Workflow 
The main methods and software are presented in the figure below:  
![alt text](https://github.com/shkot190819/serratia/blob/main/workflow_cheme.png?raw=true)

#### Uploading genomic assemblies
Get ids for all assemblies from ncbi_assembly by 'serratia' query:  
```
Entrez.email = "suhininaev@mail.ru"  
handle = Entrez.esearch(db="assembly", term='serratia', retmax='2600')  
record = Entrez.read(handle)  
ids = record['IdList']  
```
2449 ids (Date of access: 13.03.23)  

A python script [download_assemblies.py](download_assemblies.py) was written to download genomic assemblies from GenBank.

#### Reconstruction of a phylogenetic tree based on a distance matrix  
A distance matrix between genomic assemblies was constructed using [sourmash](https://github.com/sourmash-bio/sourmash).  
Sourmash is a Python package that quickly compares potentially very large sets of DNA and protein sequences. This functionality can be used to, for example, cluster transcriptomes or genomes, to identify the taxonomy of new isolate or metagenome-assembled genomes, or to determine the taxonomic composition of a new metagenome sequence by comparing it against a database of reference genomes.  
Sourmash signatures contain one or multiple sub-sampled representations of DNA or protein sequences (FracMinHash sketches). Each FracMinHash sketch contains hashes (and optionally, hash abundances) that represent a subset of k-mers from the original sequences. The sourmash sketch command consistently subsamples k-mers across different sequences, so we can compare (e.g. intersect) sketches to understand sequence similarity. The command line function sourmash compare estimates similarity and containment metrics.
```
> sourmash sketch dna -p scaled=1000,k=31 ../assemblies/*.fna.gz
> sourmash compare -p 32 ./signatures/*.sig -o ./sourmash_results/serra_cmp --distance-matrix --ksize 31 --csv ./sourmash_results/dist_matrix.csv
```

The phylogenetic tree was reconstructed using a [hierarchical clustering](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/hclust) algorithm from the r stats package and [ggtree](https://guangchuangyu.github.io/software/ggtree/) package
```
dd <- as.dist(d)     # d - distance matrix
hclust_avg <- hclust(dd, method='average')
p1 <- ggtree(hclust_avg, layout="circular") 
plot(p1)
```

##### Phylogenetic tree annotation
The resulting phylogenetic tree was annotated according to the following parameters: attributed species, nearest species, taxonomy check status, type strain. We used [ANI_report_prokaryotes.txt](https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/) from the ncbi database to obtain this data.
ANI_report_prokaryotes.txt provides Average Nucleotide Identity (ANI) data that can be used to evaluate the taxonomic identity of genome assemblies of interest.  Also included is the ANI status which GenBank uses as a basis for decisions about taxonomic identity of public genome assemblies. Specific methods used can be found [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6978984/).
This file contains ANI data for all latest archaeal & bacterial genome assemblies in GenBank, together with the paired RefSeq genome assemblies when they exist.  
Accessing the data 5 April 2023.  
Reconstruction of an annotated tree. The annotation layers were superimposed manually in Adobe Photoshop CS5, as there are too many genomes to build a clear heatmap.
```
p1 <- ggtree(hclust_avg, layout="circular") %<+% tax_table +
  geom_tippoint(aes(color=declared_species), size=6)+ 
  scale_color_manual(values = custom.col)
plot(p1)
```
The information on the type strains was taken from the [NCBI taxonomy page](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Undef&id=613&lvl=3&keep=1&srchmode=1&unlock)
After clicking on the name of the species (species, not strain), we got a link to Genome in the table at the top right of the page. Clicking on it provided the Genbank ID of the assembly for the reference strain. It is important to note that the type strains are not available for all species. We have not considered unclassified Serratia, only the major species.

![alt text](https://github.com/shkot190819/serratia/blob/main/tree.jpg?raw=true)


#### Reconstruction of a phylogenetic tree based on a core gene alignment
##### Format conversion
Annotation errors, fragmented assemblies and contamination represent a major challenge for pangenome analysis. We have used [Panaroo](https://gtonkinhill.github.io/panaroo/#/) to tackle these challenges using a sophisticated framework for error correction that leverages information across strains through a population graph-based pangenome representation.  
To run panaroo, we converted the assemblies downloaded from genbank twice. The first conversion was done with the [bp_genbank2gff3.pl script](https://manpages.ubuntu.com/manpages/focal/en/man1/bp_genbank2gff3.1p.html). This script uses Bio::SeqFeature::Tools::Unflattener and Bio::Tools::GFF to convert GenBank flatfiles to GFF3 with gene containment hierarchies mapped for optimal display in gbrowse.
```
> bp_genbank2gff3 *.gbff -o ../2gff_annotations
```
For the second conversion we used the convert_refseq_to_prokka_gff.py script, which can be found among the panaroo software scripts.  
The converted .gff files smaller than 1 Mb were deleted and not included in the further comparison. 

##### Run panaroo
Run panaroo on the converted files:
```
panaroo -i *.gff -o ./results/ --clean-mode strict -a core -t 32
```
  --clean-mode {strict,moderate,sensitive}
                        strict:
                        Requires fairly strong evidence (present in  at least
                        5% of genomes) to keep likely contaminant genes. Will
                        remove genes that are refound more often than they were
                        called originally.
  -a {core,pan}, --alignment {core,pan}  
                        Output alignments of core genes or all genes. Options
                        are 'core' and 'pan'
  --aligner {prank,clustal,mafft}
                        Specify an aligner. Options:'prank', 'clustal', and
                        default: 'mafft'

##### Clear the alignment of ambiguously defined bases 
Next, clear the alignment of ambiguously defined bases (N, U, etc., anything that is not a normal base). We replaced these bases with the most common normal base in that position throughout the alignment. If the part of the alignment has N everywhere and gaps in all other sequences, then such sections were removed. An appropriate script was provided to me to clear the alignment. Since I didn't write it myself, I didn't attach it to the repository. 
```
python3 clean_from_amb_sites.py -f core_gene_alignment_copy.aln
```
##### Clear alignment from recombination events
Recombination events can introduce incongruence between different regions of the genome, leading to misleading results in phylogenetic analyses. Therefore, it is important to identify and remove such events to obtain accurate and reliable bacterial phylogenies. Methods for detecting recombination are applied to the DNA alignment to identify regions that may have been subject to recombination, and then these regions are either removed from the analysis or treated separately to account for their potential confounding effects.  
We used [ClonalFrameML](https://github.com/xavierdidelot/ClonalFrameML), which uses maximum likelihood inference to simultaneously detect recombination in bacterial genomes and account for it in phylogenetic reconstruction.  
There are two input files needed to run ClonalFrameML. The first one is a starting tree, which must be in Newick format. The second one is an alignment of the sequences which can be either in fasta format or extended multi fasta (XMFA) format.  
The basic command for running ClonalFrameML is as follows:
```
ClonalFrameML newick_file seq_file output_prefix [OPTIONS]
```
The starting tree was generated from the alignment using [FastTree](https://bioconda.github.io/recipes/fasttree/README.html)

```
fasttree -gtr -nt core_gene_alignment_copy.fixed.fasta > fast_tree.tre
```
-gtr -- generalized time-reversible model (nucleotide alignments only)  
-nt -- fasttree  supports  fasta  or  phylip  interleaved  alignments  
By default fasttree expects protein alignments,  use -nt for nucleotides fasttree reads standard  input if no alignment file is given
                        
Polytomy was observed in the constructed tree. ClonalFrameML does not work with trees containing polytomy. To overcome this obstacle, we used the [ete3](http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html) package for python. When a tree contains a polytomy (a node with more than 2 children), the method resolve_polytomy() can be used to convert the node into a randomly bifurcated structure in which branch lengths are set to 0. This is really not a solution for the polytomy but it allows to export the tree as a strictly bifurcated newick structure, which is a requirement for some external software.

Run ClonalFrameML:
```
ClonalFrameML ../binary_tree_ete.nw ../core_gene_alignment_copy.fixed.fasta serra
```
serra - output_prefix (running ClonalFrameML produces several output files, each of which starts with the output_prefix specified in the command line)  
Output: 
Wrote inferred importation status to serra.importation_status.txt               
Wrote processed tree to serra.labelled_tree.newick                              
All done in 5308.67 minutes. 

File `importation_status.txt` contains the list of reconstructed recombination events. There is one line for each event, the first column indicates the branch on which the event was found, and the second and third columns indicate the first and last genomic positions affected by the recombination event. 

![alt text](https://github.com/shkot190819/serratia/blob/main/ClonalFrameML_output_file_.png?raw=true)

All nodes were given the same name, sorted by start and end columns, and then used the bedtools to merge all recombination ranges
```
bedtools merge -i recomb.bed > recomb_merged.txt
```
All recombinations combined into 479 ranges.  
[Bio.AlignIO](https://biopython.org/docs/1.76/api/Bio.AlignIO.html) was used to remove recombination events from the alignment:  
```
alignment = AlignIO.read(open("../core_gene_alignment_copy.fixed.fasta"), "fasta")
for record in alignment:
    for i in range(len(recomb_merged)-1,-1,-1):
        s = recomb_merged.loc[:,1][i]-1              
        e = recomb_merged.loc[:,2][i]                 
        record.seq = record.seq[:s] + record.seq[e:]
AlignIO.write(alignment, '../clonal/align_wo_rec.fasta', 'fasta')
```
The code with comments is provided in the [delete_recombinations.ipynb](delete_recombinations.ipynb) 

According to the metrics calculated with ClonalFramemML, the effect size ratio for recombination events versus mutations is 1.2. However, a very large proportion of the original alignments have been removed. 
The overestimation of regions associated with recombination is probably due to the fact that we are working with different species within a genus. Regions close to individual species are marked as recombinant, hence this method is better suited to more clonal bacterial populations. However, we have found that the level of recombination in the genus Serratia is quite low, so we do not need to remove these regions. The next part of the analysis must be repeated for alignment without removing recombinations to make a final conclusion. This step is now in the process of selecting the optimal evolutionary model. The alignment is very large, so the processing at each stage takes days. 

##### Core genome SNPs
[Snp-sites](https://github.com/sanger-pathogens/snp-sites) finds snp sites from a multi fasta alignment file. SNPs are a common type of genetic variation that occur when a single nucleotide differs between individuals or populations. SNP markers are derived from a common nucleotide alignment, which means that they are identified by comparing the DNA sequences of different individuals or populations at specific locations in the genome. Once identified, these SNPs can be used to build a phylogenetic tree that represents the evolutionary relationships among different organisms or populations. The key advantage of using SNPs to build a phylogenetic tree is that they are abundant and distributed throughout the genome, allowing for a fine-scale resolution of differences among organisms. Moreover, because SNPs can be rapidly and accurately identified using modern DNA sequencing techniques, they have become an increasingly popular tool for studying the genetic diversity and evolution of a wide range of organisms.
```
snp-sites -mvp -o align_SNP.aln align_wo_rec.fasta
```
 -m     output a multi fasta alignment file (default)  
 -v     output a VCF file  
 -p     output a phylip file  
 -o STR specify an output filename [STDOUT]  


##### Evolutionary model
Determining the optimal evolutionary model is important because it helps ensure that the resulting phylogenetic tree accurately reflects the evolutionary history of the organisms in question. Different evolutionary models make different assumptions about the rate and pattern of molecular evolution, and these assumptions can significantly affect the accuracy and reliability of the tree. Therefore, by selecting the most appropriate evolutionary model, researchers can account for differences in evolutionary rates and patterns among different branches of the tree. 
[ModelTest-NG](https://github.com/ddarriba/modeltest) is a tool for selecting the best-fit model of evolution for DNA and protein alignments. 
ModelTest-NG v0.1.7
```
modeltest-ng -i align_SNP.aln.snp_sites.aln -p 30 
```
Summary:
Partition 1/1:
| Model   | Score       | Weigth |
|---------|-------------|--------|
| TIM3+G4 | 119795.1551 | 0.9191 |
| GTR+G4  | 97567.3222  | 0.8785 |
| JC      | 27656138.0  | 1.0000 |

The JC phylogenetic model is a simple model of DNA evolution which uses a single rate parameter to describe the probability of nucleotide substitution. It assumes that all nucleotide positions are equally likely to undergo substitution and that the rate of substitution is constant throughout the evolutionary history of the sequences being analyzed. Despite its simplicity, the JC model has been widely used in phylogenetic analyses of DNA sequence data and can provide a useful baseline for more complex evolutionary models.

##### Phylogenetic tree building
[RAxML-NG](https://github.com/amkozlov/raxml-ng) is a phylogenetic tree inference tool which uses maximum-likelihood (ML) optimality criterion. Its search heuristic is based on iteratively performing a series of Subtree Pruning and Regrafting (SPR) moves, which allows to quickly navigate to the best-known ML tree.
```
raxml-ng --threads 30 --msa align_SNP.aln.snp_sites.aln --model JC
```

##### Result tree
![alt text](https://github.com/shkot190819/serratia/blob/main/20230502_tree_.jpg?raw=true)

### Conclusions
- *S. marcescens* assemblies, which has not been taxonomically verified by the NCBI, is taxonomically reassigned to the *Serratia nevei* isolated in 2020.
- *Serratia bockelmannii* is not a separate species and can be considered as *S. marcescens*

The *S. marcescens* group was found to contain a significant clade, which has been identified as *S. nevei*, suggesting a potential novel taxonomic entity within the genus *Serratia*. These observations underscore the need for a comprehensive revision of the current taxonomy of this genus.


### Literature
Cho, G. S., Stein, M., Brinks, E., Rathje, J., Lee, W., Suh, S. H., & Franz, C. M. A. P. (2020). *Serratia nevei sp. nov.* and *Serratia bockelmannii sp. nov.*, isolated from fresh produce in Germany and reclassification of *Serratia marcescens* subsp. *sakuensis* Ajithkumar et al. 2003 as a later heterotypic synonym of Serratia marcescens subsp. marcescens. Systematic and applied microbiology, 43(2), 126055. https://doi.org/10.1016/j.syapm.2020.126055  

Brown et al, (2016), sourmash: a library for MinHash sketching of DNA, Journal of Open Source Software, 1(5), 27, doi:10.21105/joss.00027  

Tonkin-Hill, G., MacAlasdair, N., Ruis, C. et al. Producing polished prokaryotic pangenomes with the Panaroo pipeline. Genome Biol 21, 180 (2020). https://doi.org/10.1186/s13059-020-02090-4  

Didelot, X., & Wilson, D. J. (2015). ClonalFrameML: efficient inference of recombination in whole bacterial genomes. PLoS computational biology, 11(2), e1004041. https://doi.org/10.1371/journal.pcbi.1004041  

Morgan N. Price and others, FastTree: Computing Large Minimum Evolution Trees with Profiles instead of a Distance Matrix, Molecular Biology and Evolution, Volume 26, Issue 7, July 2009, Pages 1641–1650, https://doi.org/10.1093/molbev/msp077  

Page AJ, Taylor B, Delaney AJ, Soares J, Seemann T, Keane JA, Harris SR. SNP-sites: rapid efficient extraction of SNPs from multi-FASTA alignments. Microb Genom. 2016 Apr 29;2(4):e000056. doi: 10.1099/mgen.0.000056  

Diego Darriba and others, ModelTest-NG: A New and Scalable Tool for the Selection of DNA and Protein Evolutionary Models, Molecular Biology and Evolution, Volume 37, Issue 1, January 2020, Pages 291–294, https://doi.org/10.1093/molbev/msz189  

Alexey M Kozlov and others, RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference, Bioinformatics, Volume 35, Issue 21, 1 November 2019, Pages 4453–4455, https://doi.org/10.1093/bioinformatics/btz305  






