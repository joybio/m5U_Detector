#!/root/miniconda3/bin/python

"""Create a summary that shows the pass through ratio"""

__date__ = "2020-9-8"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "PKU.jia.group"

#imports
import re
import os
import optparse
#sort  package
from optparse import OptionParser

parser = OptionParser('Usage: %prog ')
parser.add_option('-i','--input',
                dest='input',
                help='IGV_bedgraph')

parser.add_option('-o','--out',
                dest='out',
                help='out annotation file')

(options,args) = parser.parse_args()


data = open(options.input,"r")
out = open(options.out,"w")
out.write("tRNA_chr\ttRNA_start\ttRNA_stop\ttRNA_aa\ttRNA_codon\ttRNA_strand\ttRNA_seq\tchr\tstart\tstop\tnucletide\tcov_wt1\tcov_wt2\tcov_mut1\tcov_mut2\tpass_through_wt1\tpass_through_wt2\tpass_through_mut1\tpass_through_mut2\tpos\tFC\n")
bg_wt1 = 1
bg_wt2 = 1
bg_mut1 = 1
bg_mut2 = 1
start_pos = 1
n = 1
for i in data:
	line = i.strip()
	i = i.strip().split("\t")
	start = int(i[8])
	if start - start_pos == 1:
		n+=1
		ratio_wt1 = round((int(i[11]) + 1)/bg_wt1 * 100,1)
		ratio_wt2 = round((int(i[12]) +1)/bg_wt2 * 100,1)
		ratio_mut1 = round((int(i[13]) + 1)/bg_mut1 * 100,1)
		ratio_mut2 = round((int(i[14]) + 1)/bg_mut2 * 100,1)
		FC = round((ratio_mut1 + ratio_mut2)/(ratio_wt1 + ratio_wt2),2)
		bg_wt1 = int(i[11]) + 1
		bg_wt2 = int(i[12]) + 1
		bg_mut1 = int(i[13]) + 1
		bg_mut2 = int(i[14]) + 1
		start_pos = int(i[8])
		out.write(line + "\t" + str(ratio_wt1) + "\t" + str(ratio_wt2) + "\t" + str(ratio_mut1) + "\t" + str(ratio_mut2) + "\t" + str(n) + "\t" + str(FC) + "\n")
	elif start_pos -start == 1:
		n+=1
		ratio_wt1 = round((int(i[11]) + 1)/bg_wt1 * 100,1)
		ratio_wt2 = round((int(i[12]) +1)/bg_wt2 * 100,1)
		ratio_mut1 = round((int(i[13]) + 1)/bg_mut1 * 100,1)
		ratio_mut2 = round((int(i[14]) + 1)/bg_mut2 * 100,1)
		FC = round((ratio_mut1 + ratio_mut2)/(ratio_wt1 + ratio_wt2),2)
		bg_wt1 = int(i[11]) + 1
		bg_wt2 = int(i[12]) + 1
		bg_mut1 = int(i[13]) + 1
		bg_mut2 = int(i[14]) + 1
		start_pos = int(i[8])
		out.write(line + "\t" + str(ratio_wt1) + "\t" + str(ratio_wt2) + "\t" + str(ratio_mut1) + "\t" + str(ratio_mut2) + "\t" + str(n) + "\t" + str(FC) + "\n")
	else:
		start_pos = int(i[8])
		bg_wt1 = int(i[11]) + 1
		bg_wt2 = int(i[12]) + 1
		bg_mut1 = int(i[13]) + 1
		bg_mut2 = int(i[14]) + 1
		n = 1
		out.write(line + "\t" + "100" + "\t" + "100" + "\t" + "100" + "\t" + "100" + "\t" + str(n) + "\t" + "1" + "\n")
data.close()
out.close()





