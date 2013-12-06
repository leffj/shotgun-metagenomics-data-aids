#!/usr/bin/env python

"""
script to make gene ko codes from img-qiime package not redundant
"""

input_fp = '/Users/leffj/Dropbox/genomes_traits_project/fresh_start/KO_dictionary_img/img-qiime-25oct2012/gene_ko_pathway.txt'
output_fp = '/Users/leffj/Dropbox/genomes_traits_project/fresh_start/KO_dictionary_img/gene_ko_dictionary.txt'
outFail_fp = '/Users/leffj/Dropbox/genomes_traits_project/fresh_start/KO_dictionary_img/gene_ko_dictionary_fails.txt'

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