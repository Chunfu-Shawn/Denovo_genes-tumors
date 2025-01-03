import pandas as pd
import numpy as np
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
parser.add_argument("--sample_name", type=str, help="(Required) Sample name",
                    default="human_brain_ribo")
parser.add_argument("--txt_path", type=str, help="(Required) Output RPF ditribution file of ribotish quality (default: ribotish_qual.txt)")
parser.add_argument("--offset_path", type=str, help="(Required) Output offset file of ribotish quality (default: ribotish.para.py)")
parser.add_argument("--RPF_start_distr_file", type=str, help="(Required) Output: Distribution of RPF near start codon (considering offset)")
parser.add_argument("--RPF_stop_distr_file", type=str, help="(Required) Output: Distribution of RPF near stop codon (considering offset)")
parser.add_argument("--frame_distr_file", type=str, help="(Required) Output: Distribution of RPF 5' end across three reading frames in all annotated codons")
args = parser.parse_args()

def shift_list(ls: list, offset):
    """
    :param ls:
    :param offset:
    :return:
    """
    if offset != "dropped":
        mod = len(ls)
        ids = [[(item[0]+offset)%mod, item[1]] for item in enumerate(ls)]
        ids.sort(key=lambda item: item[0])
        return [item[1] for item in ids]
    else:
        # if no offset, output 0
        return [0] * len(ls)

# read offset file
exec(open(args.offset_path).read())

# drop the reads without offset
for l in range(27,33):
    if l not in offdict.keys():
        offdict[l] = "dropped"
# remove 'm0'
del offdict['m0']
print(offdict)

# read quality data
## txt_path:  row 1-5 (5' end match RPFs)
### row1: RPF length distribution
### row2: Distribution of RPF near start codon
### row3: Distribution of RPF near stop codon
### row4: Reads count of RPF 5' end across three reading frames in all annotated codons
### row5: The RPF profile throughout the CDS regions in 3 frame
## txt_path:  row 6-10 (5' end mismatch RPFs)

start_distr_pd = pd.DataFrame(index=range(-40,20))
stop_distr_pd = pd.DataFrame(index=range(-40,20))
frame_distr_pd = pd.DataFrame(index=range(0,3))

# merge reads with offset across different length
with open(args.txt_path, "r") as file:
    for num, line in enumerate(file):
        if num == 1:
            start_distr_dic=eval(line.strip())
            start_distr = np.sum([shift_list(start_distr_dic[i],offdict[i]) for i in offdict.keys()],
                                 axis=0, dtype=np.int32).tolist()
            start_distr_pd[args.sample_name] = start_distr

        if num == 2:
            stop_distr_dic=eval(line.strip())
            stop_distr = np.sum([shift_list(stop_distr_dic[i],offdict[i]) for i in offdict.keys()],
                                axis=0, dtype=np.int32).tolist()
            stop_distr_pd[args.sample_name] = stop_distr

        if num == 3:
            frame_distr_dic=eval(line.strip())
            frame_len = len(frame_distr_dic[25])
            frame_distr = [0] * frame_len
            frame_distr = np.sum([shift_list(frame_distr_dic[i],offdict[i]) for i in offdict.keys()],
                                 axis=0, dtype=np.int32).tolist()
            frame_distr_pd[args.sample_name] = frame_distr

# save files
start_distr_pd.T.to_csv(args.RPF_start_distr_file, index=True, sep="\t")
stop_distr_pd.T.to_csv(args.RPF_stop_distr_file, index=True, sep="\t")
frame_distr_pd.T.to_csv(args.frame_distr_file, index=True, sep="\t")
