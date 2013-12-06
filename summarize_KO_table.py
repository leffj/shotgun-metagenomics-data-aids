#!/usr/bin/env python


import sys
import os
import argparse


"""
Script to build a summary table at a specified KO hierarchy level given
a gene abundance table and a KO/corresponding hierarchy dictionary. This
script takes into account multiple hierarchies for a single KO number
by counting the abundance of the gene in all categories to which it
belongs.
"""

parser = argparse.ArgumentParser()
requiredArgs = parser.add_argument_group('required arguments')
requiredArgs.add_argument('-i', '--geneTable_fp', action='store', dest='geneTable_fp',
	required=True, help='input filepath, tab delimited')
requiredArgs.add_argument('-o', '--output_dir', action='store', dest='output_dir',
	required=True, help='output directory for summaries')
requiredArgs.add_argument('-d', '--ko_dictionary_fp', action='store', dest='ko_dictionary_fp',
	required=True, help='input filepath for KO dictionary file')
args = parser.parse_args()

geneTable_fp = args.geneTable_fp
output_dir = args.output_dir
koDictionary_fp = args.ko_dictionary_fp



def import_dictionary(dictionary_fp):
	rawDict = open(koDictionary_fp,'U')
	outDict = {}
	for line in rawDict:
		key = line.split('\t')[0]
		entry = '\t'.join(line.strip().split('\t')[1:])
		outDict[key] = entry
	rawDict.close()
	return outDict

def get_cats_from_hierarchies(geneTable,koDictionary,level):
	# get sampleIDs
	sampleIDs = geneTable.readline().strip().split('\t')[1:]
	# loop through each KO
	summaryCatsAbunds = {}
	for line in geneTable:
		splitline = line.strip().split('\t')
		ko = splitline[0]
		abunds = splitline[1:]
		abunds = [float(x) for x in abunds]
		# skip kos with abundance 0 for every sample
		if sum(abunds) == 0:
			continue
		try:
			koHiers = koDictionary[ko].split('\t')
		except KeyError:
			print 'No match in dictionary for:%s' %(ko)
			continue
		# loop through hierarchies for ko and pull cats at desired level
		# only pull the same cat once to avoid double counting
		cats = []
		for hier in koHiers:
			hierAsList = hier.strip().split(';')
			if not hierAsList[level-1].strip() in cats:
				cats.append(hierAsList[level-1].strip())
		for cat in cats:
			if not cat in summaryCatsAbunds:
				summaryCatsAbunds[cat] = abunds
			else:
				summaryCatsAbunds[cat] = [(x + y) for x, y in \
				 zip(summaryCatsAbunds[cat], abunds)]
	return sampleIDs,summaryCatsAbunds

def get_relative_abundances(summary):
	sums = []
	for key in summary:
		if sums == []:
			sums = summary[key]
		else:
			sums = [(x + y) for x, y in zip(sums,summary[key])]
	relAbunds = summary
	for key in summary:
		relAbunds[key] = [(x / y) for x, y in zip(summary[key], sums)]
	return relAbunds

def main():
	# import KO dictionary
	koDictionary = import_dictionary(koDictionary_fp)

	# create output directory if doesn't exist and create
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)

	# get output filepath prefix
	input_file_basename, input_file_ext = \
	 os.path.splitext(os.path.split(geneTable_fp)[1])

	# loop through each KO, get category for each hierarchy to which the
	# KO belongs/given level, abundances, and return summary
	levels = [1,2,3]
	for level in levels:
		sampleIDs,summary = \
		 get_cats_from_hierarchies(open(geneTable_fp,'U'),koDictionary,level)

		# get relative abundances
		summaryRAs = get_relative_abundances(summary)

		out_fp = '%s/%s_L%s.txt' %(output_dir,input_file_basename,level)
		out = open(out_fp,'w')

		out.write('Category\t' + '\t'.join(sampleIDs) + '\n')
		for key in summaryRAs:
			RAs = [str(x) for x in summaryRAs[key]]
			outLine = '%s\t%s\n' %(key,'\t'.join(RAs))
			out.write(outLine)
		out.close()


if __name__ == "__main__":
	main()