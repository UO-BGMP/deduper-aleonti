#!/usr/bin/python3

import argparse
import fileinput
import gzip
import re
import os

###################################################################################################################
######################################## DEDUPER SCRIPT - Bi624 ###################################################
###################################################################################################################



#################### Setting up argparse variables ################################################################

# Sets up required and optional argparse fields

parser = argparse.ArgumentParser(description="Python script to deduplicate single-end and paired-end sequencing reads. Requires a SAM file for input. The script will output a deduplicated SAM file, with all PCR duplicates removed. Must input a pre-sorted SAM file, which can be done using samtools 'sort' on the command line.")
parser.add_argument('-f','--file', required=True, help='the absolute file path of the SAM file for deduplication', type=str)
parser.add_argument('-p', '--paired', action='store_true', required=False, help='a flag indicating that data is paired-end; if not set, the script will assume the data is single-read')
parser.add_argument('-u','--umifile', required=False, help='an optional file containing a complete list of acceptable/expected UMI sequences', type=str)
args = parser.parse_args()

file = args.file			
pe_flag = args.paired
umifile = args.umifile
filename = os.path.basename(file)

# Prints messages to terminal depending on which flags were set by the user

if args.umifile:
	print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n NOTE: Because an UMI file was provided, UMIs must exactly match a known sequence or the read will be tossed... \n")
else:
	print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n NOTE: UMI file not provided... UMIs will be assumed to be randomers... \n")

# This script currently can't handle PE data, so if the '-p' flag is set the program will exit

if pe_flag == True:
	print(" NOTE: Input SAM file is paired-end... \n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
	raise parser.error("\n\n !!! Currently no paired-end functionality available, exiting... !!!\n")
else:
	print(" NOTE: '-p' flag not used, Input file is assumed to be single-read... \n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

################### Setting up dictionaries/lists/tuples, defining functions ########################################

umilist = []											# list to store acceptable UMIs if they are provided by the user
umicounts = {}											# counts the occurence of each UMI
readinfo = {}											# dictionary of tuples that store pertinent information about each read for deduplication

if umifile is not None:
	with open(umifile, "r") as uf:						# opens UMI file to retrieve UMI sequences (if one is provided)
		for line in uf:
			line = line.strip("\n").strip("\r")
			umilist.append(line)						# adds UMI to the list of acceptable UMIs

	for umi in umilist:
		if umi not in umicounts:						# incrimentor to keep track of how many times each UMI appears
			umicounts[umi] = 0							
		else:
			umicounts[umi] += 1

############ Function to check for the presence of an UMI in the user-provided UMI list

def umi_checker(header):
	header = header.split(":")
	if header[7] in umicounts:							# header index 7 is the UMI; if this UMI is in the user-provided UMI list, it is saved to the seq_umi variable
		seq_umi = header[7]
		return(seq_umi)
	else:
		pass											# if it isn't in the list, the line is skipped and will not be written to the output file
		
############ Function to account for soft-clipping in the start position of a read 

def position_corrector(cigar, pos):
	clipping_regex = r"[0-9]+[S]"						# this regex searches for CIGAR strings containing an "S" and grabs the numbers preceding it
	if "S" in cigar:
		softclip = str(re.search(clipping_regex, cigar).group())
		clipped_value = softclip.split("S")[0]
		correct_pos = pos - int(clipped_value)			# this takes the soft-clipped value and corrects it to report the true sequence start position
	else:
		correct_pos = pos								# if no "S"" is detected, the starting position is assumed to be correct
	return correct_pos

############ Function to check the strand of read 

def strand_determinator(bw_flag):
	if ((bw_flag & 16) != 16):			# this is obviously way oversimplified, just testing it out for now
		strand = "reverse"
	else:
		strand = "forward"
	
# 	if ((bw_flag & 1) != 1):
	# if bw_flag == 145:
	# 	strand = "reverse"
	# 	print("PE strand is reverse")   # this is only activated by -p argparse flag
	return strand	

################## Iterating through SAM file ################################################################

i = 0 													# counter variable to keep track of line number
prevChrom = 1											# variable to store the previous chromosome number; defaults to 1 at the beginning of the script

output = open(filename + "_deduped","w") 				# opens file that non-duplicate reads will be written to

with open(file, "r") as sf:
	for line in sf:
		i += 1 
		if i % 1000000 == 1:
			print("Deduping...")
			
		if line.startswith("@"):
			output.write(line)						
			pass										# skips header lines (for analysis), writes to output file
			
		else:
			clean_line = line.strip("\n").strip("\r").split("\t")
			header = clean_line[0]
			randomer = header.split(":")[7]
			bw_flag = int(clean_line[1])
			chrom = clean_line[2]
			pos = int(clean_line[3])
			cigar = clean_line[5]
			
			if umifile is not None:									# these functions will only be executed if the user provides an UMI list via the '-u' argument
				seq_umi = umi_checker(header)						# checks to make sure the UMI is free of sequencing errors
				strand = strand_determinator(bw_flag)				# Checks which strand the read is from and whether or not it's mapped
				correct_pos = position_corrector(cigar, pos)		# Ensures that any soft-clipping is accounted for prior to saving read information
				info_tuple = (seq_umi, strand, correct_pos)			# Saves all read information as a tuple
			
				if info_tuple not in readinfo:						# prints the first occurence of a line to the output file
					readinfo[info_tuple] = seq_umi
					umicounts[umi] += 1
					output.write(line)
			
				elif info_tuple in readinfo:			# skips duplicate reads based on the 'readinfo' tuple dictionary
					pass

				if prevChrom != chrom:					# resets the readinfo dictionary when a new chromosome is encountered
					readinfo = {}
					prevChrom = chrom
				
			if umifile is None:											# if an UMI list is not provided, the UMIs are assumed to be randomers
				if "N" in randomer:
					pass
					
				else: 
					strand = strand_determinator(bw_flag)				# Checks which strand the read is from and whether or not it's mapped
					correct_pos = position_corrector(cigar, pos)		# Ensures that any soft-clipping is accounted for prior to saving read information
					info_tuple = (randomer, strand, correct_pos)		# Saves all read information as a tuple
					
				if randomer not in readinfo:			# prints the first occurence of a line to the output file
					readinfo[info_tuple] = randomer
					output.write(line)
					
				elif info_tuple in readinfo:			# skips duplicate reads based on the 'readinfo' tuple dictionary
					pass

				if prevChrom != chrom:					# resets the readinfo dictionary when a new chromosome is encountered
					readinfo = {}
					prevChrom = chrom

output.close()	
print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n All done deduping! Please enjoy your shiny new deduplicated file \n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")				
					


# sample command for the terminal: python leonti_deduper.py -u STL96.txt -f /projects/bgmp/shared/deduper/Dataset1.sam
