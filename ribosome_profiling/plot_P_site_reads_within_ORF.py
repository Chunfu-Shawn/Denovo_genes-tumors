import pandas as pd
import os
import pybedtools
import argparse
import matplotlib.pyplot as plt

__author__ = "Chunfu xiao"
__contributor__="..."
__copyright__ = ""
__credits__ = []
__license__ = ""
__version__="1.0.0"
__maintainer__ = "Chunfu xiao"
__email__ = "chunfushawn@gmail.com"

# Arguments
parser = argparse.ArgumentParser(description="Arguments for this script")
parser.add_argument("--gene_id", default="gene_id", type=str, help="(Required) Gene ID", required=True)
parser.add_argument("--p_site_read_count", default="./test/test.cds.psite.depth.txt",
                    type=str, help="(Required) Read count file of P-site reads within ORF per gene", required=True)
parser.add_argument("--ORF_bed_file", default="./test/test.cds.bed",
                    type=str, help="(Required) BED file of target ORFs", required=True)
parser.add_argument("--output_dir", default="./test",
                    type=str, help="(Required) output directory", required=True)
parser.add_argument("--study", default=None,
                    type=str, nargs="+", help="Study name for filtering")
parser.add_argument("--run", default=None,
                    type=str, nargs="+", help="Sample run ID for filtering")
parser.add_argument("--start_end", default=None,
                    type=int, nargs="+", help="Start and end position of P-site reads for filtering")
args = parser.parse_args()

def count_total_P_site_reads(p_site_read_count, ORF_bed, study=None, run=None):
    for region in ORF_bed:
        chrom = region.chrom
        start = region.start
        end = region.end
        strand = region.strand
        p_site_read_count_position = {pos: 0 for pos in range(int(end)-int(start))}
        # filter by positions and samples
        p_site_read_count_region = p_site_read_count[(p_site_read_count.Chrom == chrom) & \
                                                     (p_site_read_count.Pos > int(start)) & (p_site_read_count.Pos <= int(end))]
        if study != None:
            p_site_read_count_region = p_site_read_count_region[p_site_read_count_region.Study.isin(study)]
        if run != None:
            p_site_read_count_region = p_site_read_count_region[p_site_read_count_region.Run.isin(run)]
        print(p_site_read_count_region)
        for i in range(len(p_site_read_count_region)):
            if strand == "+":
                position = p_site_read_count_region.iloc[i, 1] - 1 - start
            else:
                position = end - p_site_read_count_region.iloc[i, 1] 
            # read counts of relative positions
            p_site_read_count_position[position] += p_site_read_count_region.iloc[i, 2]
    return p_site_read_count_position

def plot_P_site_reads(p_site_read_count_position, output_file, log_y=False):
    # plot the read counts within ORF using different colors of frame 0, 1, 2
    plt.figure(figsize=(9, 4))
    plt.bar(p_site_read_count_position.keys(), p_site_read_count_position.values(), \
            width=1, color=["#08306B","#BABABA","#E0E0E0"])
    plt.xlabel('Position of ORF', fontsize=11,)
    plt.ylabel('# of P-site reads', fontsize=11)
    if log_y:
        plt.yscale("log")
    plt.subplots_adjust(left=0.1, bottom=0.15)
    # save pdf
    plt.savefig(output_file)
    
def plot_three_nucleotide_periodicity(p_site_read_count_position, output_file, start_end=None):
    read_counts_frame={pos: 0 for pos in range(3)}
    total_reads=1
    # count reads in frame 0, 1, 2
    for pos in p_site_read_count_position.keys():
        # if filter reads by start and end position within ORF
        if isinstance(start_end, list) and len(start_end)==2:
            total_reads += p_site_read_count_position[pos] if pos>=start_end[0] and pos<=start_end[1] else 0
            read_counts_frame[pos%3] += p_site_read_count_position[pos] if pos>=start_end[0] and pos<=start_end[1] else 0
        else:
            total_reads += p_site_read_count_position[pos]
            read_counts_frame[pos%3] += p_site_read_count_position[pos]
    plt.figure(figsize=(3, 6))
    plt.bar(read_counts_frame.keys(), [x/total_reads for x in read_counts_frame.values()], \
            width=0.8, color=["#08306B","#BABABA","#E0E0E0"])
    plt.xlabel('Reading frame', fontsize=11,)
    plt.ylabel('Proportion of P-site reads', fontsize=11)
    plt.margins(0.1)
    plt.subplots_adjust(left=0.25, bottom=0.15)
    # save pdf
    plt.savefig(output_file)

def main():
    # read data
    p_site_read_count = pd.read_csv(
        args.p_site_read_count, 
        sep="\t", header=0)
    ORF_bed = pybedtools.BedTool(args.ORF_bed_file)
    # summrized read counts per position
    p_site_read_count_position = count_total_P_site_reads(p_site_read_count, ORF_bed, args.study, args.run)
    # plot 
    plot_P_site_reads(
        p_site_read_count_position, 
        args.output_dir + "/number_of_P-site_reads_within_ORF."+ args.gene_id + ".pdf")
    plot_three_nucleotide_periodicity(
        p_site_read_count_position, 
        args.output_dir + "/three_nt_periodicity_within_ORF."+ args.gene_id + ".pdf",
        args.start_end)

if __name__=="__main__":
    print(args.study, args.run, args.start_end)
    main()