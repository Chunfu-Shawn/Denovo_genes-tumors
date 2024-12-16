################################################
#File Name: run.sh
#Author: rbase    
#Mail: xiaochunfu@stu.pku.edu.cn
#Created Time: Wed 06 Dec 2023 09:13:45 PM CST
################################################

#!/bin/sh 

usage(){
  echo "Usage: bash $(basename $0) -1 label -2 workDir -3 geneAA -4 querydb -5 msData"
  echo "Author: rbase"
  echo "Description: This script perform peptide searching for MS spectrum by blastp and verify the unique matches."
  echo "Date: 2024-7-8"
  echo "------------------------------------------------"
  echo "OPTIONS"
  echo -e "     --work_dir [Mandatory] \t\t\t\tThe work and output path"
  echo -e "     --orf_bed_file [Mandatory] \t\t\t\tBED file with ORFs"
  echo -e "     --maf_dir [Mandatory] \t\t\t\tFolder with MAF alignments, separated by chrm (chrm1.maf, chrm2.maf...)"
  echo -e "     --output_dir [Mandatory] \t\t\t\tOutput folder name"
  echo -e "     --nwk_file [Mandatory] \t\t\t\tnwk alignment"
  echo -e "     --pep_fasta_file [Mandatory] \t\t\t\tFASTA with ORF proteins"
  echo -e "     -h|--help\t\t\t\tprint this help page"
  exit
}

####################################################
## De novo status evaluation based on 120 mammals ##
####################################################

mkdir -p $work_dir/tmp/results_120
mkdir -p $work_dir/results_120
        
echo "Calculating multiple alignments. Recommend to parallelize the searches for multiple ORFs."
python3 1_extract_multiple_alignments.py \
    -b $orf_bed_file -m $maf_dir -o $work_dir/results_120 -f yes

echo "Calculating ancestral sequences and estimnate intact ORF ancestrally"
bash 2.1_ancestral_sequences.sh $work_dir/results_120 \
    $nwk_file \
    $orf_bed_file.120mammals \
    $pep_fasta_file

echo "Performing similarity searches across orthologous regions to trace protein age"
python3 3_sequence_specificity.py --prot_dir $work_dir/results_120 --prot_tar $pep_fasta_file

echo "Curating these de novo ORFs manually using primate synteny data"
### maltiple sequence alignments for these ORFs and define "common disablers"