# NGS-Indel_Coder
NGS-Indel_Coder Package including 11 python scripts, manual and example files.

NGS-Indel Coder v.1.0.0 Manual<br/>
By Julien Boutte, PhD<br/>
Shannon Straub Lab (Hobart and William Smith Colleges, Geneva, NY, USA)<br/>
June 2019<br/>

NGS-Indel Coder was developed to detect and omit false positive indels inferred from assemblies of short read sequence data. This tool, divided in five parts proposed several options. This pipeline used several tools including 2MATRIX (Salinas and Little 2014) to code indels as binary characters and BLAST to detect exon positions (Altschul et al. 1990). Output files were generated for IQ-TREE software (Nguyen et al. 2015, Chernomor et al. 2016). Nevertheless, NGS-Indel Coder output fasta files can be used with any software coding indels using aligned fasta files. <br/>

List of python scripts:<br/>

1-parsing_Samtools_depth-files.py<br/>
2-Indel_validation.py<br/>
3-Indel_validation.py<br/>
4-Indel_validation.py<br/>
5-indel_deletion.py<br/>
6-IQTREE_binary_matrices_creation.py<br/>
7-IQTREE_DNA_matrices_creation.py<br/>
8-nexus_files_creation.py<br/>
9-identification_boundaries.py<br/>
10-partitioned_nexus_files_creation.py<br/>
11-Delete_small_partitions.py<br/>

When using NGS-Indel Coder please cite:<br/>

NGS-Indel Coder: A pipeline to code indel characters in phylogenomic data with an example of its application in milkweeds (Asclepias), Julien Boutte, Mark Fishbein, Aaron Liston, and Shannon C.K. Straub. In press in MPE.<br/>

The NGS-Indel_Coder_Manual.pdf file contains:<br/>

I.	About NGS-Indel Coder, citation<br/>
II.	Downloading NGS-Indel Coder, getting help<br/>
III.	Input files format<br/>
IV.	NGS-Indel Coder command lines and options<br/>
V.	Output files<br/>
VI.	Example<br/>
VII.	References<br/>
