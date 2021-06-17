#!/root/miniconda3/bin/python
__date__ = "2021-5-20"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "PKU.jia.group"


import re
import os
import optparse
from optparse import OptionParser


parser = OptionParser("Usage: %prog -i tRNA.merge.depth.bed -o tRNA.strand.merge.depth.bed")
parser.add_option("-i","--input",dest = "input",
		help = "input file: tRNA.merge.depth.bed") 
parser.add_option("-o","--output",dest = "out",
		help = "Output file: tRNA.strand.merge.depth.bed.")

(options,args) = parser.parse_args()

out = open(options.out,"w")

with open(options.input,"r") as f:
	for i in f:
		i = i.strip().split("\t")
		out.write(i[6] + "\t" + i[7] + "\t" + i[8] + "\t" + i[3] + "\t" + i[4] + "\t" + i[5] + "\n")
out.close()



