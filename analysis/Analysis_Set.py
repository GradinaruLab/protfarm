import Sequence_Library

class Analysis_Set:

	def __init__(self):

		self.sequence_libraries = {}

	def add_library(self, file, name):

		sequence_library = Sequence_Library(file)
		self.sequence_libraries[name] = sequence_library


	# libraries_of_interest: list of identifiers for the library that you want the specificity of
	# libraries_to_compare: list of identifiers for the library that you want to compare the specificity against
	#
	# Returns: Nx2 matrix, 1st column is sequence, 2nd column is specificity score
	def get_specificity(self, libraries_of_interest_name, libraries_to_compare_names):

		specificity_dict={}
		library_of_interest = self.sequence_libraries[libraries_of_interest_name]
		num_comparing_libraries = len(libraries_to_compare_names)
		for key in library_of_interest:
			comparing_libraries_sum = 0
			library_of_interest = library_of_interest.get_sequence_counts()
			for libraries_compare_index in len(libraries_to_compare_names):
				compare_library_name = libraries_to_compare_names[libraries_compare_index]
				comparing_library = self.sequence_libraries[compare_library_name]
				comparing_library = comparing_library.get_sequence_counts()
				comparing_libraries_sum += comparing_library[key]/sum(comparing_library.values())
			library_of_interest_counts = library_of_interest[key]/sum(library_of_interest.values())
			comparing_libraries_sum += libraries_of_interest_counts
			sequence_specificity = library_of_interest_counts/(comparing_libraries_sum*num_comparing_libraries)
			specificity_dict[key] = sequence_specificity

		for libraries_compare_index in len(libraries_to_compare_names):
			compare_library_name = libraries_to_compare_names[libraries_compare_index]
			comparing_library = self.sequence_libraries[compare_library_name]
			comparing_library = comparing_library.get_sequence_counts()
			for key in comparing_library:
				if key not in specificity_dict:
					specificity_dict[key] = 0




	def get_enrichment(self, libraries_of_interest, libraries_to_compare):

		return 0