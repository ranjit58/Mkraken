#! /bin/python

# This script will merge reports from the kraken pipeline into a single report and can filter and sort results as .OTU.

import sys
import pandas as pd
import numpy as np
from collections import OrderedDict as OD
from os import listdir, makedirs
from os.path import isfile, isdir, join, exists
import argparse

__author__="Brody DeSilva and Ranjit Kumar"
__email__="bdesilva@uab.edu, rkumar@uab.edu"

# change the behavior to 
	# output a folder with the merged file and each of the filtered files
		# check for existence of the folder in the first place, if exist do not overwrite
	# have a -input_suffix clause to specify files (ex: .report)
	# 


# specify input and output files and optionally, filter id
parser = argparse.ArgumentParser(description="""Merge Kraken files into a single file with an additional average column. 
	Save the file as (input)_merge_kraken.txt Optionally, filter this file based on a classification. Save the sorted file
	(input)_merge_kraken_(classification character).OTU into an OTU table format with the average column in descending 
	order of percent abundance.""")
parser.add_argument('-s','--suffix',help="""Provide the suffix argument to select files within the input folder
provided.""")
parser.add_argument('-f','--filter',help="""Filter the merged_kraken file using any of the following characters to 
specify the taxonomic classification [D K P C O F G S U - no]. D : Domain, K : Kingdom, P : Phylum, C : Class, O : Order, 
F : Family, G : Genus, S : Species, U : Unclassified, - : anything that is not listed, no : do not filter.""", choices=['D','K','P','C','O','F','G','S','U','-', 'no'])
parser.add_argument('-a', '--avg_type',help="""Use this to specify using the Kraken average  precision of 2 decimal places
.""",action='store_true')
parser.add_argument('input',help="""Specify a folder with at least 2 files in it; make sure only the desired files are in
the folder.""")
parser.add_argument('output',help="""Specify an output file prefix. This prefix will be asserted at the beginning of the 
any file printed from this folder. Required output: (input)_merge_kraken.txt. Optional Output: 
(input)_merge_kraken_(classification_character).OTU""")
parser.add_argument('file_format',help="""Specify a file with the complete taxonomy listing of the Kraken Database. This file is
used to sort the output files to have the original kraken hierarchy of taxonomy.""")
parser.add_argument('-v','--verbose', help="""Use this argument to receive additional output from the program.""",
action='store_true')
args = parser.parse_args()

if args.verbose:
	import time
	timein = time.time()

files = []
# check if input is file or folder, if folder get all files in a list
if isdir(args.input):
	files = [''.join([args.input,'/',f]) for f in listdir(args.input) if isfile(join(args.input, f))] 
	if len(files) < 2:
		sys.exit('At least 2 files are needed in a folder.\nAborting.\n') # quit execution
	if args.verbose:
		print('There were %d files found.\n') % len(files)
else:
	sys.exit('A folder was not provided, or the path given was incorrect. Check and try again.\n') # quit execution

# apply the suffix to the list of files
if args.suffix != None:
	files = [f for f in files if f[-len(suffix):] == suffix]
	if len(files) < 1:
		sys.exit("""The suffix provided did not apply to any files within the folder given. Please check the suffix and 
		input folder and try again.\n""") # quit execution

# check for and create output folder
if exists(args.output):
	sys.exit("""The folder specified already exists, please input a new folder name or remove the current folder and
	try again.\n""") # quit execution
else:
	makedirs(args.output)	

# check for the file_format file
if isfile(args.file_format):
	pass
else:
	sys.exit("""Check the path to the file format file and try again. No file was found with the name/path provided""")

# specify the average type to be provided in the output files
if args.avg_type:
	REALNUMS = 0
else:
	REALNUMS = 1


# import tab-delimited table files
# columns of interest are the 0th, 3rd, 4th and 5th
# create a new column with the file name for each new file
master_list = []
temp_list = []
flag = 1
file_names = []
ffn = args.file_format # format file name
files.append(ffn)

def removeFileEndings(f): # Recursively remove file endings from the files so  tables display sample names clearly
	while len(f.split('.')) != 1 or len(f.split('/')) != 1:
		f=f.split('.')[0]
		f=f.split('/')[-1]		
	return(f)
if args.verbose:
	print('Reading in files and analysing...\n')


for file in files: #  for each file
	if file != ffn:
		file_names.append(removeFileEndings(file)) #remove file endings
	f = pd.read_csv(file, sep='\t', header=None) #  open file as .csv using the panda module
	TOTAL = float(f.iloc[0][1] + f.iloc[1][1])/100 # UNSURE IF THIS IS ALWAYS TRUE
	# parse through each file and find key-value pairs based on the taxa id and the taxa
	if REALNUMS != 1:
		f.columns = ['per', 'b', 'c', 'classif', 'tax_id', 'tax']
		a = f.loc[:,['per', 'classif', 'tax_id', 'tax']] # get percent abundance, classification, tax_id, and taxonomy
	else:
		f.columns = ['b', 'per', 'c', 'classif', 'tax_id', 'tax']
		a = f.loc[:,['per', 'c', 'classif', 'tax_id', 'tax']] # get percent abundance, classification, tax_id, and taxonomy
	if flag == 1:  # if this is the first run, then create the master series to be compared and added
		# create dict with 3 keys
		if REALNUMS != 1:
			master_list = dict(zip(zip(a['tax_id'].values, a['tax'].values, a['classif'].values), a['per'].values)) 
		else:
			master_list = dict(zip(zip(a['tax_id'].values, a['tax'].values, a['classif'].values),
			a['per'].values/TOTAL))
		m_df = pd.Series(master_list)  # create the master series with index value a tuple of the 3 keys
		if args.verbose:
			print('(rows, columns) -> (elements, files)')
			print(m_df.shape)
		flag = 0 # ensure continuing on to temporary files
	else: # do the exact same thing with the temporary series (i.e.temporary because each file after the first one)
		# create dictionary
		if REALNUMS != 1:
			temp_list = dict(zip(zip(a['tax_id'].values, a['tax'].values, a['classif'].values), a['per'].values)) 
		else:
			temp_list = dict(zip(zip(a['tax_id'].values, a['tax'].values, a['classif'].values),
			a['per'].values/TOTAL)) 
		temp_series = pd.Series(temp_list) # create series from dict
		if file == ffn:
			if REALNUMS != 1:
				temp_list = OD(zip(zip(a['tax_id'].values, a['tax'].values, a['classif'].values), a['per'].values))
			else:
				temp_list = OD(zip(zip(a['tax_id'].values, a['tax'].values, a['classif'].values),
				a['per'].values/TOTAL))
			temp_series = pd.Series(temp_list) # create series from dict
			break
		m_df = pd.concat([m_df, temp_series], axis=1) # add line to the dataframe
		if args.verbose:
			print(m_df.shape) # show number of unique keys (as rows) and the number of files read (as columns)
file_names.append('avg')
m_df = pd.concat([m_df, m_df.mean(axis=1)], axis=1) # add element average across all the files as a column
					 	    # if a sample does not have an element present, it will not be used for the average
m_df.set_axis(1, file_names) # set the columns to be the names of the files (done by order of reading)
m_df_sort = m_df.sort_values('avg',0,False) # sort
m_df_tax = m_df.reindex(temp_series.index)
m_df_tax = m_df_tax.dropna(0,'all')

fid = join(args.output, 'mkr_sort_abd.txt') # file output name : Merged Kraken Report, the one sorted by abundance
if args.verbose:
	print('Writing file %s...\n') % fid
with open(fid, 'w') as f: # open file with write permissions (overwrite any data previous)
	m_df_sort.to_csv(f, '\t') # write a tab-delimited file using pandas' built-in method

fid = join(args.output, 'mkr_taxonomy.txt') # file output name : Merged Kraken Report, the one with the set taxonomy
if args.verbose:
	print('Writing file %s...\n') % fid
with open(fid, 'w') as f: # open file with write permissions (overwrite any data previous)
	m_df_tax.to_csv(f, '\t') # write a tab-delimited file using pandas' built-in method


def filterAndSort(m_df, filterspec, output, verbose):
	tax_df = m_df.loc[m_df.index[:].get_level_values(2) == filterspec] # filter out the matching taxonomy
	# put in the format of the standardized OTU table
	# taxid :  samples :  taxonomy
	multi = tax_df.index
	tax_df = tax_df.set_index(multi.get_level_values(0))
	tax_df['taxonomy'] = multi.get_level_values(1)
	tax_df['taxonomy'] = [x.strip() for x in tax_df['taxonomy'].values] # remove spaces at beginning of taxonomy names
	tax_df = tax_df.sort_values('avg',0,False) # sort the dataframe based on the magnitude of the avg descending

	fid = join(output, ''.join(['mkr_', filterspec,'.OTU']))
	if verbose:
		print('Writing OTU table to output file %s...\n') % fid
	with open(fid, 'w') as f: # open file with write permissions (overwrite any data previous)
			tax_df.to_csv(f, '\t') # write a tab-delimited file using pandas' built-in method

# if filtering
if args.filter != 'no':
	if args.verbose:
		print('Filter provided: %s') % args.filter
		print('Sorting the merged file...\n')
	if args.filter == None:
		for ident in ['D','K','P','C','O','F','G','S','U','-']:
			filterAndSort(m_df, ident, args.output, args.verbose)
	else:
			filterAndSort(m_df, ident, args.output, args.verbose)

	if args.verbose:
		timeout = time.time()
		print('Time Elapsed : %d') % (timeout-timein)
		print(timeout-timein)
else:
	if args.verbose:
		print('No filter provided. Program finished successfully')
