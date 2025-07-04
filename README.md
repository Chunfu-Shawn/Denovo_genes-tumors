# Part of analysis code in project: Oncogenic Roles of Young Human De Novo Genes and Their Potential as Neoantigens in Cancer Immunotherapy
## evolution_orf: Analyzing the evolution of ORF sequences and structures across primate and mammalian species
1. 1_extract_multiple_alignments: extract multiple alignments of ORF sequences
2. 2.1_ancestral_sequences.sh, 2.2_parsing_ancestors.py: reconstruct ancestral sequences
3. 3_sequence_specificity.py: evaluate the sequence conservation of ORFs and trace similar proteins across the included species. 
These scripts were adjusted from https://github.com/jorruior/riboseq_orfs_analyses.

## ribosome_profiling: Analyzing the 3-nt periodicity around de novo ORFs
1. parse_ribotish_qual.py: extract and merge read information from results generated by RiboTISH tool
2. generate_P_site_reads.py: generate bam file containing peptidyl site (P-site) reads, using offset file generated from RiboTISH tool
3. plot_P_site_reads_within_ORF.py: plot read counts of P-site reads for each position of ORF, using the depth file (format according to the test file)