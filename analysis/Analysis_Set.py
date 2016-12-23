from analysis.Sequence_Library import Sequence_Library
import csv
from utils import DNA
import math
from workspace import Workspace as ws
from workspace import Database as db

class Analysis_Set:

	def __init__(self):

		self.sequence_libraries = {}
		self.sequence_length = 0

	def add_library(self, library):
		if isinstance(library,str):
			library = db.get_library(library)

		sequence_library = Sequence_Library(library)
		self.sequence_length = sequence_library.get_sequence_length()
		self.sequence_libraries[library.name] = sequence_library

	def get_sequence_length(self, by_amino_acid = True):

		if by_amino_acid:
			return self.sequence_length / 3
		else:
			return self.sequence_length

	def get_libraries(self):
		return self.sequence_libraries

	# libraries_of_interest: list of identifiers for the library that you want the specificity of
	# libraries_to_compare: list of identifiers for the library that you want to compare the specificity against
	#
	# Returns: Nx2 matrix, 1st column is sequence, 2nd column is specificity score
	def get_specificity(self, library_of_interest_names, libraries_to_compare_names, by_amino_acid = True, count_threshold = 10):

		specificity_dict={}

		library_of_interest_total_count = 0;

		if isinstance(library_of_interest_names, list):

			sequence_counts = {}

			for library_of_interest_name in library_of_interest_names:
				library_of_interest = self.sequence_libraries[library_of_interest_name]
				library_of_interest_total_count += library_of_interest.get_total_count()
				library_of_interest_counts = library_of_interest.get_sequence_counts(by_amino_acid, count_threshold = 0)

				for sequence, count in library_of_interest_counts.items():
					if sequence not in sequence_counts:
						sequence_counts[sequence] = count
					else:
						sequence_counts[sequence] += count

			below_threshold_sequences = set()
			for sequence, count in sequence_counts.items():
				if count < count_threshold:
					below_threshold_sequences.add(sequence)

			for sequence in below_threshold_sequences:
				del(sequence_counts[sequence])

			num_comparing_libraries = len(libraries_to_compare_names) + len(library_of_interest_names)
		else:
			library_of_interest = self.sequence_libraries[library_of_interest_names]
			library_of_interest_total_count = library_of_interest.get_total_count()
			sequence_counts = library_of_interest.get_sequence_counts(by_amino_acid, count_threshold = count_threshold)

			num_comparing_libraries = len(libraries_to_compare_names) + 1

		libraries_to_compare = []
		libraries_to_compare_total_counts = []

		for library_to_compare_name in libraries_to_compare_names:
			library_to_compare = self.sequence_libraries[library_to_compare_name]
			libraries_to_compare_total_counts.append(library_to_compare.get_total_count())
			library_to_compare = library_to_compare.get_sequence_counts(by_amino_acid, count_threshold)
			libraries_to_compare.append(library_to_compare)

		#print "Comparing to " + str(num_comparing_libraries) + " library(ies)"
		#print "Library of interest total count: " + str(library_of_interest_total_count)

		for key in sequence_counts:

			comparing_libraries_presence = 0

			for library_to_compare_index in range(len(libraries_to_compare)):
				library_to_compare = libraries_to_compare[library_to_compare_index]

				library_presence = 0
				if key in library_to_compare:
					#print key + " is in! Count: " + str(library_to_compare[key])
					library_presence = library_to_compare[key] * 1.0 /libraries_to_compare_total_counts[library_to_compare_index]
					comparing_libraries_presence += library_presence
				else:
					pass
					#print key + " is NOT in!"

				#print "Presence of " + key + " in " + libraries_to_compare_names[library_to_compare_index] + ": " + str(library_presence)

			library_of_interest_presence = 1.0*sequence_counts[key]/library_of_interest_total_count
			comparing_libraries_presence += library_of_interest_presence

			sequence_specificity = library_of_interest_presence/(comparing_libraries_presence/num_comparing_libraries)
			specificity_dict[key] = sequence_specificity

		for libraries_compare_index in range(len(libraries_to_compare_names)):
			compare_library_name = libraries_to_compare_names[libraries_compare_index]
			comparing_library = self.sequence_libraries[compare_library_name]
			comparing_library = comparing_library.get_sequence_counts(by_amino_acid = by_amino_acid, count_threshold = count_threshold)
			for key in comparing_library:
				if key not in specificity_dict:
					specificity_dict[key] = 0

		return specificity_dict

	def get_enrichment(self, library_of_interest_names, starting_library_names, by_amino_acid = True, count_threshold = 10, Log_Scale=True, zero_count_magic_number = 0.9,include_zero_count=False):

		if isinstance(library_of_interest_names, list):

			library_of_interest_total_count = 0
			sequence_counts = {}

			for library_of_interest_name in library_of_interest_names:
				library_of_interest = self.sequence_libraries[library_of_interest_name]
				library_of_interest_total_count += library_of_interest.get_total_count()
				library_of_interest_counts = library_of_interest.get_sequence_counts(by_amino_acid, count_threshold = 0)

				for sequence, count in library_of_interest_counts.items():
					if sequence not in sequence_counts:
						sequence_counts[sequence] = count
					else:
						sequence_counts[sequence] += count

			below_threshold_sequences = set()
			for sequence, count in sequence_counts.items():
				if count < count_threshold:
					below_threshold_sequences.add(sequence)

			for sequence in below_threshold_sequences:
				del(sequence_counts[sequence])

			library_of_interest = sequence_counts
		else:
			library_of_interest = self.sequence_libraries[library_of_interest_names]
			library_of_interest_total_count = library_of_interest.get_total_count()
			library_of_interest = library_of_interest.get_sequence_counts(by_amino_acid, count_threshold = 0)

		if isinstance(starting_library_names, list):

			aggregate_starting_library = {}

			for starting_library_name in starting_library_names:
				starting_library = self.sequence_libraries[starting_library_name]
				starting_library_total_count = starting_library.get_total_count()
				starting_library = starting_library.get_sequence_counts(by_amino_acid, count_threshold = count_threshold)

				for sequence, sequence_count in starting_library.items():

					if sequence not in aggregate_starting_library:
						aggregate_starting_library[sequence] = sequence_count
					else:
						aggregate_starting_library[sequence] += sequence_count

			starting_library = aggregate_starting_library

		else:
			starting_library = self.sequence_libraries[starting_library_names]
			starting_library_total_count = starting_library.get_total_count()
			starting_library = starting_library.get_sequence_counts(by_amino_acid, count_threshold = count_threshold)

		enrichment_dict = {}
		for sequence in starting_library:

			if sequence not in library_of_interest or library_of_interest[sequence] == 0:
				library_of_interest[sequence] = zero_count_magic_number

			if (Log_Scale):
				fold_enrichment = (library_of_interest[sequence]* 1.0 / library_of_interest_total_count) / (starting_library[sequence] * 1.0/ starting_library_total_count)
				enrichment_dict[sequence] = math.log10(fold_enrichment)
			else:
				fold_enrichment = (library_of_interest[sequence]* 1.0 / library_of_interest_total_count) / (starting_library[sequence] * 1.0/ starting_library_total_count)
				enrichment_dict[sequence] = fold_enrichment

		# For the sequences that don't exist in the starting library
		if (include_zero_count):
			for sequence in library_of_interest:
				if sequence in starting_library:
					continue

				fold_enrichment = (library_of_interest[sequence] * 1.0 / library_of_interest_total_count) / (zero_count_magic_number / starting_library_total_count)

				if Log_Scale:
					enrichment_dict[sequence] = math.log10(fold_enrichment)
				else:
					enrichment_dict[sequence] = fold_enrichment

		return enrichment_dict

	def get_enrichment_library_of_interest(self, library_of_interest_name, starting_library_name, by_amino_acid = True, count_threshold = 10, Log_Scale=False, zero_count_magic_number = 0.9):
		library_of_interest = self.sequence_libraries[library_of_interest_name]
		library_of_interest_total_count = library_of_interest.get_total_count()
		library_of_interest = library_of_interest.get_sequence_counts(by_amino_acid, count_threshold = 0)
		starting_library = self.sequence_libraries[starting_library_name]
		starting_library_total_count = starting_library.get_total_count()
		starting_library = starting_library.get_sequence_counts(by_amino_acid, count_threshold = count_threshold)
		enrichment_dict = {}
		for sequence in library_of_interest:
			if sequence not in starting_library:
				fold_enrichment = (library_of_interest[sequence] * 1.0 / library_of_interest_total_count) / (zero_count_magic_number / starting_library_total_count)

			if Log_Scale:
				enrichment_dict[sequence] = math.log10(fold_enrichment)
			else:
				enrichment_dict[sequence] = fold_enrichment
		return enrichment_dict

	def export_enrichment(self, filename, starting_libary_name, \
		by_amino_acid = False, count_threshold = 0, log_scale = False):

		num_sequence_libraries = len(self.sequence_libraries)

		cumulative_counts = {}
		cumulative_enrichments = {}

		library_index = 0

		header_row = ['Sequence']
		if not by_amino_acid:
			header_row.append('Amino Acid')

		for library_name, library in self.sequence_libraries.items():

			print('Getting counts and enrichment for ' + library_name)

			header_row.append(library_name)
			header_row.append(library_name + ' enrichment')

			if isinstance(starting_libary_name, dict):
				fold_enrichments = self.get_enrichment(library_name, starting_libary_name[library_name], by_amino_acid, count_threshold = 0, Log_Scale = log_scale, include_zero_count = True)
			else:
				fold_enrichments = self.get_enrichment(library_name, starting_libary_name, by_amino_acid, count_threshold = 0, Log_Scale = log_scale, include_zero_count = True)
			
			library_counts = library.get_sequence_counts(by_amino_acid, count_threshold = 0, filter_invalid = False)

			for sequence, sequence_count in library_counts.items():

				if sequence not in cumulative_counts:
					cumulative_counts[sequence] = [0] * num_sequence_libraries
					cumulative_enrichments[sequence] = [0] * num_sequence_libraries

				if sequence in fold_enrichments:
					cumulative_enrichments[sequence][library_index] = fold_enrichments[sequence]

				cumulative_counts[sequence][library_index] = sequence_count

			library_index += 1

		data = []

		for sequence, sequence_counts in sorted(cumulative_counts.items(), key=lambda kv: sum(kv[1]), reverse=True):

			sequence_row = [sequence]

			if not by_amino_acid:
				sequence_row.append(DNA.translate_dna_single(sequence))

			if sum(sequence_counts) < count_threshold:
				continue

			library_index = 0

			fold_enrichments = cumulative_enrichments[sequence]

			for sequence_count in sequence_counts:
				sequence_row.append(sequence_count)
				sequence_row.append(fold_enrichments[library_index])
				library_index += 1

			data.append(sequence_row)

		ws.export_csv(filename, header_row, data)

	def export_enrichment_specificity(self, filename, starting_libary_name, \
		libraries_to_compare_names, by_amino_acid = False, \
		count_threshold = 0):

		num_sequence_libraries = len(self.sequence_libraries)

		cumulative_counts = {}
		cumulative_enrichments = {}
		cumulative_specificities = {}

		library_index = 0

		header_row = ['Sequence']
		if not by_amino_acid:
			header_row.append('Amino Acid')

		for library_name, library in self.sequence_libraries.items():

			header_row.append(library_name)
			header_row.append(library_name + ' enrichment')
			header_row.append(library_name + ' specificity')

			fold_enrichments = self.get_enrichment(library_name, starting_libary_name, by_amino_acid, count_threshold = 0, Log_Scale = False)
			specificities = self.get_specificity(library_name, libraries_to_compare_names, by_amino_acid, count_threshold = 0)
			library_counts = library.get_sequence_counts(by_amino_acid, count_threshold = 0, filter_invalid = False)

			for sequence, sequence_count in library_counts.items():

				if sequence not in cumulative_counts:
					cumulative_counts[sequence] = [0] * num_sequence_libraries
					cumulative_enrichments[sequence] = [0] * num_sequence_libraries
					cumulative_specificities[sequence] = [0] * num_sequence_libraries

				if sequence in fold_enrichments:
					cumulative_enrichments[sequence][library_index] = fold_enrichments[sequence]

				if sequence in specificities:
					cumulative_specificities[sequence][library_index] = specificities[sequence]

				cumulative_counts[sequence][library_index] = sequence_count

			library_index += 1

		data = []

		for sequence, sequence_counts in sorted(cumulative_counts.items(), key=lambda kv: sum(kv[1]), reverse=True):

			sequence_row = [sequence]

			if not by_amino_acid:
				sequence_row.append(DNA.translate_dna_single(sequence))

			if sum(sequence_counts) < count_threshold:
				continue

			library_index = 0

			fold_enrichments = cumulative_enrichments[sequence]
			specificities = cumulative_specificities[sequence]

			for sequence_count in sequence_counts:
				sequence_row.append(sequence_count)
				sequence_row.append(fold_enrichments[library_index])
				sequence_row.append(specificities[library_index])
				library_index += 1

			data.append(sequence_row)

		ws.export_csv(filename, header_row, data)
