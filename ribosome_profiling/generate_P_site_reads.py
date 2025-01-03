import pandas as pd
import os
import uuid
import array
import pybedtools
import pysam
import argparse

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
parser.add_argument("--bam_file", type=str, help="(Required) BAM file of Ribo-seq reads mapped into genome")
parser.add_argument("--ORF_bed_file", type=str, help="(Required) BED file of target ORFs")
parser.add_argument("--offset_path", type=str, help="(Required) Output offset file of ribotish quality (default: ribotish.para.py)")
parser.add_argument("--p_site_bam_file", type=str, help="(Required) Output: site-corrected reads using offset files")
args = parser.parse_args()

def shift_reads(bam, ORF_bed, offdict, tmp_path, BAM=True):
    tmp_file = pysam.AlignmentFile(tmp_path, "wb" if BAM else "w", header=bam.header)
    # drop reads without offset
    # consider positive or negative strand
    for read in bam:
        # check read length
        if read.query_length in offdict.keys():
            offset = offdict[read.query_length]
            # position
            read.pos = read.pos + offset if ORF_bed[0].strand == "+" else read.pos + read.query_length - 1 - offset
            # 1 nt sequence
            read.query_sequence = read.query_sequence[offset-1:offset] if ORF_bed[0].strand == "+" else read.query_sequence[-(offset+1):-offset]
            if read.query_qualities is None:
                read.query_qualities = array.array('B', [38])
            else:
                read.query_qualities = read.query_qualities[offset-1:offset] if ORF_bed[0].strand == "+" else read.query_qualities[-(offset+1):-offset]
            read.cigartuples = [(0, 1)]
            # add
            # print(read.pos, read.cigartuples, read.query_sequence, read.query_qualities, read.query_length)
            tmp_file.write(read)

def filter_save_reads_from_SAM_by_bed(tmp_path, ORF_bed, output_path):
    tmp_sam = pysam.AlignmentFile(tmp_path, "r")
    output_bam = pysam.AlignmentFile(output_path, "wb", header=tmp_sam.header)
    for read in tmp_sam:
        chrom = read.reference_name
        start = read.pos
        end = start + read.query_length
        read_bed = pybedtools.BedTool(f"{chrom}\t{start}\t{end}", from_string=True)
        # if overlapped
        if read_bed.intersect(ORF_bed, u = True):
            output_bam.write(read)
    tmp_sam.close()
    output_bam.close()

def filter_save_reads_from_BAM_by_bed(tmp_path, ORF_bed, output_path):
    tmp_bam = pysam.AlignmentFile(tmp_path, "rb")
    output_bam = pysam.AlignmentFile(output_path, "wb", header=tmp_bam.header)
    for region in ORF_bed:
        chrom = region.chrom
        start = region.start
        end = region.end
        for read in tmp_bam.fetch(chrom, start, end):
            output_bam.write(read)
    tmp_bam.close()
    output_bam.close()

def main():
    output_dir = ".".join(args.p_site_bam_file.split(".")[:-3])
    tmp_path = output_dir + '.exon.psite.sam'
    # read ribo-seq reads
    bam = pysam.AlignmentFile(args.bam_file, "rb")
    # read ORF bed
    ORF_bed = pybedtools.BedTool(args.ORF_bed_file)
    # shift reads and keep P-sites
    shift_reads(bam, ORF_bed, offdict, tmp_path, BAM=False)
    # keep reads within ORFs
    filter_save_reads_from_SAM_by_bed(tmp_path, ORF_bed, args.p_site_bam_file)
    bam.close()
    


if __name__=="__main__":
    # read offset file
    if os.path.exists(args.offset_path):
        exec(open(args.offset_path).read())
        # remove 'm0'
        del offdict['m0']
        if len(offdict.keys()) != 0:
            print(offdict)
            main()
        else:
            print("offsets are absent, pass")
    else:
        print("offset file doesn't exits, pass")