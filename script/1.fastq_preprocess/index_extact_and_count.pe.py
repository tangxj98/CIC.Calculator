#!/usr/bin/env python3 

import os
import sys
import gzip
import pprint

#debug = True
debug = False

#primer_indexes = ['CATGC', 'ACGAG', 'GTACT', 'CGATC', 'ATCGT', 'GCATG', 'TAGCA', 'CAACG', 'GCTCA', 'TAGGT', 'GAGCA', 'TGTCT', 'AGTGG', 'GTGCC', 'CTAGA', 'CGGAT', 'GACGG']


#primer_indexes = ['TGCTA', 'CATGC', 'ACGAG', 'GTACT', 'CGATC', 'ATCGT', 'GCATG', 'TAGCA', 'CAACG', 'GCTCA', 'TAGGT', 'GAGCA', 'TGTCT', 'AGTGG', 'GTGCC', 'CTAGA', 'CGGAT', 'GACGG', 'TCAAC', 'ATGAA']

primer_indexes = ['AGATCG']

#barcode_indexes = ['TCGA', 'AGCT', 'GATC', 'CTAG', 'ACTG']
barcode_indexes = ['ACTG', 'AGCT', 'CTAG', 'GATC', 'GTCA', 'TCGA']


# the start of this sequence is found 5' of the primer index approximately 25 bases away
primer_constant_5prime_alpha = "AGCAG"
primer_constant_5prime_alpha_distance = 31

# the start of this sequence is found 5' of the primer index approximately 5 bases away
primer_constant_5prime_beta = "GCTAGTGCTAG"
primer_constant_5prime_beta_distance = 11


barcode_3prime_constant = "AAGTA"


def extract_index(line=None):

	index_set = dict()

	if line is None:
		return index_set

	found_primer_indexes = dict()

	found_barcode_indexes = dict()

	found_primer_constants_alpha = find_all_positions(line=line, index=primer_constant_5prime_alpha)

	found_primer_constants_beta = find_all_positions(line=line, index=primer_constant_5prime_beta)

	for primer_index in primer_indexes:
		position_set = find_all_positions(line=line, index=primer_index)
		if position_set:
			found_primer_indexes[primer_index] = position_set


	for barcode_index in barcode_indexes:
		position_set = find_all_positions(line=line, index=barcode_index + barcode_3prime_constant)
		if position_set:
			found_barcode_indexes[barcode_index] = position_set







	#if debug:
	#	print(line + " start")

	for barcode_index in found_barcode_indexes:
		for barcode_index_position in found_barcode_indexes[barcode_index]:

			lowest_primer_index = None
			lowest_primer_index_position = len(line)
			lowest_padding_value = len(line)

			for primer_index in found_primer_indexes:
				
				for primer_index_position in found_primer_indexes[primer_index]:

					primer_constant_5prime_alpha_found = False
					primer_constant_5prime_beta_found = False
					total_padding = 0 

					for primer_constant_5prime_alpha_location in found_primer_constants_alpha:

						padding = abs(primer_index_position - primer_constant_5prime_alpha_location - primer_constant_5prime_alpha_distance)
						if padding < 3:
							primer_constant_5prime_alpha_found = True
							#if debug:
							#	print(" " * primer_constant_5prime_alpha_location + str(primer_constant_5prime_alpha) + " alpha ")
							total_padding = padding

					for primer_constant_5prime_beta_location in found_primer_constants_beta:

						padding = abs(primer_index_position - primer_constant_5prime_beta_location - primer_constant_5prime_beta_distance)

						if padding < 3:
							#if debug:
							#	print(" " * primer_constant_5prime_beta_location + str(primer_constant_5prime_beta) + " beta ")
							primer_constant_5prime_beta_found = True
							total_padding = total_padding + padding

					if not primer_constant_5prime_alpha_found or not primer_constant_5prime_beta_found:
						continue

					if total_padding < lowest_padding_value:
						lowest_padding_value = total_padding
						lowest_primer_index_position = primer_index_position
						lowest_primer_index = primer_index
			

			distance = lowest_primer_index_position - barcode_index_position

			if lowest_primer_index is not None:

				middle_index = line[barcode_index_position + len(barcode_index): barcode_index_position + len(barcode_index) + 33]
				index2 = line[ lowest_primer_index_position + 34 : lowest_primer_index_position + 44]
				if debug:
					print( line + " found" )
					print(" " * primer_constant_5prime_alpha_location + str(primer_constant_5prime_alpha) + " alpha ")
					print(" " * primer_constant_5prime_beta_location + str(primer_constant_5prime_beta) + " beta ")
					print( (" " * barcode_index_position) + barcode_index + (" " * (distance - len(barcode_index) ) ) + lowest_primer_index  + " (" + str(distance) + ")")
					print(" " * (barcode_index_position + len(barcode_index)) + middle_index )
					print(" " * (lowest_primer_index_position + 34) + index2 + " i7" )
				#return ",".join([barcode_index, middle_index, lowest_primer_index])
				return ",".join([barcode_index, middle_index, index2])


	if debug:
		print(line + "not found")
	print(line)



def find_all_positions(line=None, index=None):

	positions = list()

	if line is None or index is None:
		return positions

	search_start = 0
	no_further_instances = False

	while not no_further_instances:	

		position = None

		position = line.find(index, search_start)

		if position == -1:
			no_further_instances = True
			continue

		positions.append(position)

		search_start = position + 1

	return positions









input_fastq = sys.argv[1]

print(input_fastq)

fastq = gzip.open(input_fastq, mode='r')

count = 1 

index_set = dict()


for line in fastq:
	
	if count == 2:

		new_index = extract_index(line=line.decode().strip())

		if new_index in index_set:
			index_set[new_index] = index_set[new_index] + 1
		else:
			index_set[new_index] = 1

	count = count + 1 

	if count == 5:
		count = 1


for key in sorted(index_set, key=index_set.get, reverse=True):
  print(key, index_set[key])

