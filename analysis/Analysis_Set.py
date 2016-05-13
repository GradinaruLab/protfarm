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
	def get_specificity(self, libraries_of_interest, libraries_to_compare):

		return 0

	def get_enrichment(self, libraries_of_interest, libraries_to_compare):

		return 0