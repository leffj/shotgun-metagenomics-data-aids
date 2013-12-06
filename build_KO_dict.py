#!/usr/bin/env python

"""
script to make gene ko codes from img-qiime package not redundant
"""

import argparse

parser = argparse.ArgumentParser()
requiredArgs = parser.add_argument_group('required arguments')
requiredArgs.add_argument('-i', '--input_fp', action='store', dest='input_fp',
	required=True, help='input filepath, tab delimited')
requiredArgs.add_argument('-o', '--output_fp', action='store', dest='output_fp',
	required=True, help='output filepath for KO dictionary')
requiredArgs.add_argument('-f', '--outFail_fp', action='store', dest='outFail_fp',
	required=True, help='output filepath for KO entries that failed to be added to dictionary')
args = parser.parse_args()

input_fp = args.input_fp
output_fp = args.output_fp
outFail_fp = args.outFail_fp

input = open(input_fp,'U')

ontologyDict = {}
output = ''
outputFail = ''
for line in input:
	# if KO number hasn't been seen yet, record into dict
	hierarchies = line.strip().split('\t')
	try:
		geneID = hierarchies[1].split('; ')[3]
	except IndexError:
		outputFail = '%s\n' %hierarchies[1]
	try:
		ontologyDict[geneID]
	except KeyError:
		ontologyDict[geneID] =  '\t'.join(hierarchies[1:])

#need to fix output
for key in ontologyDict:
	output = '%s%s\t%s\n' %(output,key,ontologyDict[key])

out = open(output_fp,'w')
out.write(output)
out.close()

out2 = open(outFail_fp,'w')
out2.write(outputFail)
out2.close()