#!/usr/bin/python3

import argparse
import fileinput
import gzip

###################################################################################################################
######################################## DEDUPER SCRIPT - Bi624 ###################################################
###################################################################################################################


#################### setting up ARGPARSE variables ########################################################

# Sets up  required and optional argparse fields

parser = argparse.ArgumentParser(description="Python script to deduplicate single-end and paired-end sequencing reads. Requires a SAM file for input. The script will output a deduplicated SAM file, will all PCR duplicates removed.")
parser.add_argument('-f','--file', required=True, help='the absolute file path of the SAM file for deduplication', type=str)
parser.add_argument('-p', '--paired', action='store_true', help='a flag indicating that data is paired-end; if not set, the script will assume the data is single-read')
parser.add_argument('-u','--umifile', required=False, help='an optional file containing a complete list of acceptable/expected UMI sequences', type=str)
args = parser.parse_args()

file = args.file			
pe_flag = args.paired
umifile = args.umifile	

# Include the following argparse options
# -f, --file: required arg, absolute file path
# -p, --paired: optional arg, designates file is paired end (not single-end)
# -u, --umi: optional arg, designates file containing the list of UMIs (unset if randomers instead of UMIs)
# -h, --help: optional arg, prints a USEFUL help message (see argparse docs)
# If your script is not capable of dealing with a particular option (ex: no paired-end functionality), 
# your script should print an error message and quit


################### Set up dictionaries/lists/tuples, define functions ######################################################

umilist = []
umicounts = {}

with open(umifile, "r") as uf:						# open UMI file to retrieve UMI sequences
	for line in uf:
		line = line.strip("\n").strip("\r")
		umilist.append(line)

for umi in umilist:
	if umi not in umicounts:						# store correct index pairs in a dictionary with the value	
		umicounts[umi] = 0							# (value = count) set to 0
	else:
		umicounts[umi] += 1
		
print("This is the umi list:", umilist)

############ Function to check for the presence of an UMI in the user-provided UMI list #############



def umi_checker(umicounts):
	header = line[1]
	header = line.strip("\n").strip("\r").split(":")
	seq_umi = header[8]
	if seq_umi in umicounts:
		umicounts[seq_umi] += 1
		print(umicounts[seq_umi])
	else:
		pass
		
print(umicounts)


############ Function to account for soft-clipping in the start position of a read #############

def position_corrector(n):
	cigar = line[6]
	pos = line[4]
	if "S" in cigar:
		clipping = number associated with clipped value
		correct_pos = pos - clipping
    elif:
       correct_pos = pos
    return correct_pos
#     
# def strand_determinator(?):
# 	strand = line[2]
# #	if strand = 4 or 16, can't remember which it is:
# #		strand = reverse
# #	elif:
# #		strand = forward	
# 	return strand	
	
readinfo = {}
	
    
    
    




################## Open files and assign information to data structures ###################################



################## Iterate through SAM file ################################################################

# i = 0 
# with open(file, "r") as sf:
# 	for line in sf:
# 		i += 1 
# 		if i % 4 == 1:
# 			line = line.split("\t")
# 			header = line[0]
# 			currentChromo = [2]
# 			 if currentChromo 

# only runs umichecker function if the -u argparse flag has been set
# if umifile == "-u":		