#!/root/miniconda3/bin/python
""" if tRNA locate in the "-" strand, then reverse tRNA seq """
__date__ = "2021-5-20"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "PKU.jia.group"

import re
import os
import optparse
from optparse import OptionParser

parser = OptionParser("Usage: %prog -i tRNA.merge.seq.slice.depth.bed -o new_tRNA.merge.seq.slice.depth.bed")
parser.add_option("-i","--input",dest = "input",
		help = "input file: tRNA.merge.seq.slice.depth.bed")
parser.add_option("-o","--output",dest = "out",
		help = "Output file: new_tRNA.merge.seq.slice.depth.bed")

(options,args) = parser.parse_args()

data = open(options.input,"r")
out = open(options.out,"w")

fwd = []
rvs = []
for i in data:
	line = i
	i = i.strip().split("\t")
	if i[5] == "+":
		fwd.append(line)
	else:
		rvs.append(line)
for i in fwd:
	out.write(i)
rvs_r = rvs[::-1]
for j in rvs_r:
	out.write(j)
data.close()
out.close()
	
		





