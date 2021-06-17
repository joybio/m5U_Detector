#!/root/miniconda3/bin/python
__date__ = "2021-5-20"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "PKU.jia.group"

import re
import os
import optparse
from optparse import OptionParser


parser = OptionParser("Usage: %prog -d tRNA.merge.depth.bed -f hg38.tRNA.seq.fa -s hg38.tRNA.slice.strand.fa -o tRNA.merge.seq.depth.bed")
parser.add_option("-d","--depth",dest = "depth",
		help = "tRNA.merge.depth.bed")
parser.add_option("-f","--fa",dest = "fa",
		help = "hg38.tRNA.seq.fa")
parser.add_option("-s","--slice",dest = "slice",
		help = "hg38.tRNA.slice.strand.fa")
parser.add_option("-o","--output",dest = "out",
                help = "Output file: tRNA.merge.seq.depth.bed.")

(options,args) = parser.parse_args()

data = open(options.depth,"r")
out = open(options.out,"w")
tRNA = open(options.fa,"r")
tRNA_slice = open(options.slice,"r")

#format of tRNA (tRNA positon and seq): >chromosome:start-stop(strand)
#dict: key = chromosome:start-stop; value = nucletide sequence
tRNA_dict = {}
for i in tRNA:
	if i.startswith(">"):
		pre_key = i.lstrip(">").strip()
		key = pre_key.split("(")[0]
	else:
		tRNA_dict[key] = i.strip()
#print(tRNA_dict)
#format of tRNA_slice (single position and single nucletide): >chromosome:start-stop(strand)
#dict: key = chromosome:start-stop; value = nucletide
tRNA_slice_dict = {}
for i in tRNA_slice:
	if i.startswith(">"):
		pre_key = i.lstrip(">").strip()
		key = pre_key.split("(")[0]
	else:
		tRNA_slice_dict[key] = i.strip()
#print(tRNA_slice_dict)
#depth of single nucletide of tRNA (include tRNA sequence)
for i in data:
	i = i.strip().split("\t")
	key_tRNA = i[0]+":"+i[1]+"-"+i[2]
	key_tRNA_slice =i[6]+":"+i[7]+"-"+i[8]
	if key_tRNA in tRNA_dict.keys():
		if key_tRNA_slice in tRNA_slice_dict.keys():
			out.write(i[0] + "\t" + i[1] + "\t" +i[2] + "\t" + i[3] + "\t" +i[4] + "\t" + i[5] + "\t" + tRNA_dict[key_tRNA] + "\t" +i[6] + "\t" +i[7] + "\t" +i[8] + "\t" + tRNA_slice_dict[key_tRNA_slice] + "\t" + i[9] + "\t" +i[10] + "\t" +i[11] + "\t" +i[12] + "\n")
data.close()
out.close()
tRNA.close()
tRNA_slice.close()




