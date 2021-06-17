#!/root/miniconda3/bin/python

__date__ = "2021-5-20"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "PKU.jia.group"

import re
import os
import optparse
from optparse import OptionParser


parser = OptionParser("Usage: %prog -i homo-mature-tRNA.format.fa -o hg38.tRNA.strand.bed")
parser.add_option("-i","--input",dest = "input",
		help = "Input file; homo-mature-tRNA.format.fa. Download from: http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-mature-tRNAs.fa")
parser.add_option("-o","--output",dest = "out",
		help = "Output file; bedfile.")

(options,args) = parser.parse_args()

out=open(options.out,'w')
with open (options.input,'r') as data:
	for line in data:
		if line.startswith(">"):
			row = line.strip().split(" ")
			pos = row[12].lstrip("chr").split(":")
			print(pos)
			chrom = pos[0]
			start_stop = pos[1].replace("-","\t")
			anticodon = row[5].lstrip("(").rstrip(")")
			strand = row[13].lstrip("(").rstrip(")")
			out.write(chrom + "\t" + start_stop + "\t" + row[4] + "\t" + anticodon + "\t" + strand +'\n')
		else:
			pass
out.close()


