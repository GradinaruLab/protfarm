from utils import DNA
from workspace import Workspace as ws
from fileio import csv_wrapper

class Sequence_Library:

    # Loads the sequence library file associated with a library
    def __init__(self, library):

        alignment_file_name = ws.get_alignment_file_name(library, \
            ws.get_active_alignment())

        self._sequence_UUID_counts = csv_wrapper.read_csv_file(\
            alignment_file_name)

    def get_sequence_length(self):

        return len(self._sequence_UUID_counts[0][0])

    def get_total_count(self):
        """Returns total count of sequences"""

        sum = 0

        for sequence_UUID_count in self._sequence_UUID_counts:
            sum += sequence_UUID_count[2]

        return sum

    @property
    def sequence_UUID_counts(self):
        return self._sequence_UUID_counts

    def get_sequence_counts(self, by_amino_acid=True, count_threshold=10, filter_invalid=True):
        """Returns an Nx2 matrix, 1st column is sequence, 2nd column is count"""

        sequence_counts = {}

        for sequence_UUID_count in self._sequence_UUID_counts:

            sequence = sequence_UUID_count[0]

            if (filter_invalid or by_amino_acid) and sequence.find('N') != -1:
                continue

            sequence_count = sequence_UUID_count[2]

            if by_amino_acid:

                sequence = DNA.translate_dna_single(sequence)

                if filter_invalid and sequence.find("#") != -1:
                    continue

            if sequence not in sequence_counts:
                sequence_counts[sequence] = sequence_count
            else:
                sequence_counts[sequence] += sequence_count

        masked_sequence_counts = {}
        for sequence, sequence_count in sequence_counts.items():

            if sequence_count < count_threshold:
                pass
            else:
                masked_sequence_counts[sequence] = sequence_count

        return masked_sequence_counts

    def eliminate_bias(self,predicted_bias_percentage,num_nucleotides_off):
    """Returns an Nx2 matrix, 1st column is sequence, 2nd column is count"""

    # Initiliaze variables
    filtered_sequences = []
    current_high_count_sequence = 0;

    # Get a dictionary of unique sequences and their respective counts
    sequence_dict = self.get_sequence_counts(sequence_lib,10,True)

    #Generate an ordered list
    sorted_sequence_list = sorted(sequence_dict.items(),key=lambda x: x[1], reverse=True)


    #Iterate through the list and eliminate sequences with bias
    counter = 0
    for sequence_index in xrange(len(sorted_sequence_list)):
        if (sequence_index > len(sorted_sequence_list)):
            break
        else:       
        current_sequence_compare_value = predicted_bias_percentage*sorted_sequence_list[sequence_index][1]
        current_sequence = sorted_sequence_list[sequence_index][0]
        counter = counter + 1

        list_of_sequences_to_delete = {}
        for sequence_under_test_index in xrange(len(sorted_sequence_list)):
            if (sequence_comparable(sorted_sequence_list[sequence_under_test_index+counter][0],current_sequence,num_nucleotides_off):
                if (sorted_sequence_list[sequence_under_test_index+counter][1] <= current_sequence_compare_value):
                    # DEBUG #
                    print 'Sequence with high count is '+current_sequence
                    print 'High count sequence value is '+str(current_sequence_compare_value)
                    print 'The sequence that is the bias of the high count sequence is '+sorted_sequence_list[sequence_under_test_index+counter][0]
                    print 'That sequence value is '+str(sorted_sequence_list[sequence_under_test_index+counter][1])

                    # END DEBUG #
                    list_of_sequences_to_delete.add(sequence_under_test_index+counter)


    # Iterate in reverse through sequence list to delete bad indices
    for bad_sequence_index in list_of_sequences_to_delete[::-1]:
        sorted_sequence_list.remove(list_of_sequences_to_delete[bad_sequence_index])
            if (list_of_sequences_to_delete[bad_sequence_index] == counter):
                counter = counter - 1

    filtered_sequence_dict = dict(itertools.izip_longest(*[iter(l)] sorted_sequence_list 2, fillvalue=""))

    return filtered_sequence_dict



    def sequence_comparable(sequence_to_compare,sequence_to_be_compared_to,threshold):
    """Returns true if two string are comparable given user defined threshold"""
    sequence_letter_difference = sum(1 for x,y in zip(sequence_to_compare, sequence_to_be_compared_to) if x != y)
    if (sequence_letter_difference <= threshold):
        return True
    else:
        return False
