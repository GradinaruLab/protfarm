from utils import DNA
from workspace import Workspace as ws
from fileio import csv_wrapper
import itertools
import math

class Sequence_Library:

    # Loads the sequence library file associated with a library
    def __init__(self, library):

        alignment_file_name = ws.get_alignment_file_name(library, \
            ws.get_active_alignment())

        print("Reading CSV file for %s" % library.name)
        self._sequence_UUID_counts = csv_wrapper.read_csv_file(\
            alignment_file_name)

        self._has_UUIDs = True

        if len(self._sequence_UUID_counts) > 0:
            if math.isnan(self._sequence_UUID_counts[0][1]):
                self._has_UUIDs = False

        print("Read CSV file for %s" % library.name)



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

        if not by_amino_acid:
            sequence_counts = {x[0]: x[2] for x in self._sequence_UUID_counts}
        else:
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

        print('Eliminating bias from ' + str(len(self._sequence_UUID_counts)) + ' sequences')
        # Initialize variables
        filtered_sequences = []
        current_high_count_sequence = 0;

        # Get a dictionary of unique sequences and their respective counts
        sequence_dict = self.get_sequence_counts(False,0,True)

        #Generate an ordered list
        sorted_sequence_list = sorted(sequence_dict.items(),key=lambda x: x[1], reverse=True)


        #Iterate through the list and eliminate sequences with bias
        counter = 0
        num_filtered = 0
        duplicate_indices = 0
        list_of_sequence_indices_to_delete = set() #make it a set since you might have duplicates
        for sequence_index in range(len(sorted_sequence_list)):
            # avoids counting sequences that would yield count*percentage less than 1 (there is no sequences less than 0.99)
            if sorted_sequence_list[sequence_index][1] * predicted_bias_percentage < 1:
                break
            else:       
                current_sequence_compare_value = predicted_bias_percentage*sorted_sequence_list[sequence_index][1]
                current_sequence = sorted_sequence_list[sequence_index][0]
                counter = counter + 1
                #print('Sequence under test is '+current_sequence) 
                #print('Sequence value under test is '+str(current_sequence_compare_value))  
                
                for sequence_under_test_index in range(len(sorted_sequence_list) - counter):
                    if (Sequence_Library.sequence_comparable(sorted_sequence_list[sequence_under_test_index+counter][0],current_sequence,num_nucleotides_off)):                 
                        if (sorted_sequence_list[sequence_under_test_index+counter][1] <= current_sequence_compare_value):
                            if sequence_under_test_index+counter in list_of_sequence_indices_to_delete:
                                duplicate_indices = duplicate_indices + 1
                            # DEBUG #
                            # print('Sequence with high count is '+current_sequence)
                            # print('High count sequence value is '+str(sorted_sequence_list[sequence_index][1]))
                            # print('The sequence that is the bias of the high count sequence is '+sorted_sequence_list[sequence_under_test_index+counter][0])
                            # print('That sequence value is '+str(sorted_sequence_list[sequence_under_test_index+counter][1]))
                            num_filtered = num_filtered + 1
                            # END DEBUG #
                            list_of_sequence_indices_to_delete.add(sequence_under_test_index+counter)


        # print(str(len(list_of_sequence_indices_to_delete)))
        # print(str(duplicate_indices))        
        # Reverse sort the set 
        list_of_sequence_indices_to_delete = sorted(list_of_sequence_indices_to_delete,reverse=True)
        for bad_sequence_index in list_of_sequence_indices_to_delete:
            del sorted_sequence_list[bad_sequence_index]

        self._sequence_UUID_counts = []
        for index in range(len(sorted_sequence_list)):
            self._sequence_UUID_counts.append([sorted_sequence_list[index][0], '', sorted_sequence_list[index][1]])

        print(str(len(self._sequence_UUID_counts)) + ' sequences remaining after filtering')

    def sequence_comparable(sequence_to_compare,sequence_to_be_compared_to,threshold):
        """Returns true if two string are comparable given user defined threshold"""
        sequence_letter_difference = sum(1 for x,y in zip(sequence_to_compare, sequence_to_be_compared_to) if x != y)
        if (sequence_letter_difference <= threshold):
            return True
        else:
            return False