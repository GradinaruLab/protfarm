from Sequence_Library import Sequence_Library
import csv
import DNA as DNA

class Analysis_Set:

	def __init__(self):

		self.sequence_libraries = {}

	def add_library(self, filename, name):

		sequence_library = Sequence_Library(filename)
		self.sequence_libraries[name] = sequence_library


	# libraries_of_interest: list of identifiers for the library that you want the specificity of
	# libraries_to_compare: list of identifiers for the library that you want to compare the specificity against
	#
	# Returns: Nx2 matrix, 1st column is sequence, 2nd column is specificity score
	def get_specificity(self, library_of_interest_name, libraries_to_compare_names, by_amino_acid = True, count_threshold = 10):

		specificity_dict={}
		library_of_interest = self.sequence_libraries[library_of_interest_name]
		library_of_interest = library_of_interest.get_sequence_counts(by_amino_acid, count_threshold)

		libraries_to_compare = []
		libraries_to_compare_total_counts = []

		for library_to_compare_name in libraries_to_compare_names:
			library_to_compare = self.sequence_libraries[library_to_compare_name]
			library_to_compare = library_to_compare.get_sequence_counts(by_amino_acid, count_threshold)
			libraries_to_compare.append(library_to_compare)
			libraries_to_compare_total_counts.append(sum(library_to_compare.values()))

		library_of_interest_total_count = sum(library_of_interest.values())

		num_comparing_libraries = len(libraries_to_compare_names) + 1

		for key in library_of_interest:

			comparing_libraries_presence = 0

			for library_to_compare_index in range(len(libraries_to_compare)):
				library_to_compare = libraries_to_compare[library_to_compare_index]

				if key in library_to_compare:
					comparing_libraries_presence += library_to_compare[key]/libraries_to_compare_total_counts[library_to_compare_index]

			library_of_interest_presence = 1.0*library_of_interest[key]/library_of_interest_total_count
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

	def get_enrichment(self, library_of_interest_name, starting_library_name, by_amino_acid = True, count_threshold = 10, Log_Scale=False):

		library_of_interest = self.sequence_libraries[library_of_interest_name]
		library_of_interest_total_count = library_of_interest.get_total_count()
		library_of_interest = library_of_interest.get_sequence_counts(by_amino_acid, count_threshold = 0)
		starting_library = self.sequence_libraries[starting_library_name]
		starting_library_total_count = starting_library.get_total_count()
		starting_library = starting_library.get_sequence_counts(by_amino_acid, count_threshold = count_threshold)

		enrichment_dict = {}

		for sequence in starting_library:

			if sequence not in library_of_interest:
				enrichment_dict[sequence] = 0.0
			else:
				fold_enrichment = (library_of_interest[sequence]* 1.0 / library_of_interest_total_count) / (starting_library[sequence] * 1.0/ starting_library_total_count)
				enrichment_dict[sequence] = fold_enrichment

		return enrichment_dict

	def print_to_file(self, filename, starting_libary_name = '2_SV', by_amino_acid = False, count_threshold = 0):

		num_sequence_libraries = len(self.sequence_libraries)

		output_file = open(filename, 'w')

		writer = csv.writer(output_file)

		cumulative_counts = {}
		cumulative_enrichments = {}

		library_index = 0

		header_row = ['Sequence']
		if not by_amino_acid:
			header_row.append('Amino Acid')

		for library_name, library in self.sequence_libraries.iteritems():

			header_row.append(library_name)
			header_row.append(library_name + ' enrichment')

			fold_enrichments = self.get_enrichment(library_name, starting_libary_name, by_amino_acid, count_threshold = 0)
			library_counts = library.get_sequence_counts(by_amino_acid, count_threshold = 0, filter_invalid = False)
			
			for sequence, sequence_count in library_counts.iteritems():

				if sequence not in cumulative_counts:
					cumulative_counts[sequence] = [0] * num_sequence_libraries
					cumulative_enrichments[sequence] = [0] * num_sequence_libraries

				if sequence in fold_enrichments:
					cumulative_enrichments[sequence][library_index] = fold_enrichments[sequence]

				cumulative_counts[sequence][library_index] = sequence_count

			library_index += 1

		writer.writerow(header_row)

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

			writer.writerow(sequence_row)

		output_file.close()