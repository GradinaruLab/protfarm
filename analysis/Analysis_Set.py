from Sequence_Library import Sequence_Library

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
			comparing_library = comparing_library.get_sequence_counts()
			for key in comparing_library:
				if key not in specificity_dict:
					specificity_dict[key] = 0

		return specificity_dict

	def get_enrichment(self, libraries_of_interest, libraries_to_compare):

		return 0