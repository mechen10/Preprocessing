#!/bin/python

#######################################
# "Converting OTU Table to LeFSe compatible file"
# V2017.04.13
# Melissa Chen


#######################################
import argparse
import os
import sys

parser = argparse.ArgumentParser(
	description="Combined OTU Table with metadata for LEfsE")
parser.add_argument(
	'-i',
	'--input_OTUTable',
	help = "OTU Table from QIIME-- biom format",
	required = True,)
parser.add_argument(
	'-m',
	'--metadata',
	help = 'Mapping File with traits',
	required = True)
parser.add_argument(
	'-c',
	'--category',
	help = 'Categories (comma separated) to include in Lefse',
	type = str,
	required = True)
parser.add_argument(
	'-o',
	'--output_file',
	help = 'Name of output file',
	required = True)
	
args = parser.parse_args()

OTUFP = args.input_OTUTable
metadataFP = args.metadata
category = args.category
output_file = args.output_file

categories = category.split(",")

#######################################
# LOAD DATA

# Make OTU Table
os.system('biom convert -i ' + OTUFP + ' --to-tsv -o OTU_Table_text.txt')

# Load OTU Table
biomOpen = open("OTU_Table_text.txt", 'r') # This is in Unix format
biomTemp = []
for i in biomOpen:
	tempLine = i.strip()
	tempLineSplit = tempLine.split("\t")
	biomTemp.append(tempLineSplit)
biomOpen.close()

OTUTable = {} # Now make dictionary
for lineN in range(len(biomTemp)):
	if lineN == 0: # This is first line; skip this
		pass
	elif lineN == 1: # This is the header line
		headers = biomTemp[lineN][1:]
	else:
		OTUTable[str(biomTemp[lineN][0])] = biomTemp[lineN][1:]

# Get total counts for each site
totalCounts = {}
for h in range(len(headers)):
	totalCounts[headers[h]] = 0
for OTU in OTUTable:
	for h in range(len(headers)):
		totalCounts[headers[h]] += float(OTUTable[OTU][h])
		
os.system("rm OTU_Table_text.txt")
# 		
# Load metadata
metadataOpen = open(metadataFP, 'U') # U is for 'Universal read'-- automatically turns into Unix LF
metadataTemp = []
for i in metadataOpen:
	lineTemp = i.strip()
	lineTempSplit = lineTemp.split("\t")
	metadataTemp.append(lineTempSplit)
metadataOpen.close()

metadata = {}
metadataSites = []
for lineN in range(len(metadataTemp)):
	if lineN == 0:
		headerList = metadataTemp[lineN]
		for headerName in metadataTemp[lineN]:
			metadata[headerName] = {}
	else:
		for i in range(1,len(metadataTemp[lineN])):
			metadataSites.append(metadataTemp[lineN][0])
			sortHeader = headerList[i]
			metadata[sortHeader][metadataTemp[lineN][0]] = metadataTemp[lineN][i]
			
#######################################
# Check to make sure all names in OTU table are in metadata

if any([False for i in headers if i in metadataSites]):
	sys.exit("NOT ALL SITES IN OTU TABLE ARE IN METADATA-- PLEASE CHECK FILES AND TRY AGAIN")
	
#######################################
# Now, add metadata into file and print out

# Print headers
toWrite = "SampleID"
for i in headers:
	toWrite += "\t" + i
toWrite += "\n"

# Print data
for cat in categories:
	toWrite += str(cat)
	for i in headers:
		toWrite += "\t" + metadata[cat][i]
	toWrite += "\n"

# Print OTU abundances
for OTU in OTUTable.keys():
	toWrite += OTU
	for siteN in range(len(headers)):
		toWrite += "\t" + str(float(OTUTable[OTU][siteN])/float(totalCounts[headers[siteN]]))
	toWrite += "\n"


LEFSEfile = open(str(output_file), 'w')
LEFSEfile.write(toWrite)
LEFSEfile.close()


		

			

	





