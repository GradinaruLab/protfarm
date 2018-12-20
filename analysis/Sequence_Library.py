from utils import DNA
from workspace import Workspace as ws
from fileio import csv_wrapper
import itertools
import math
import pandas
import numpy
from utils import Sequence_Trie

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
            if isinstance(self._sequence_UUID_counts[0][1], str):
                self._has_UUIDS = True
            elif math.isnan(self._sequence_UUID_counts[0][2]):
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

    def get_sequence_counts(self, by_amino_acid=True, count_threshold=10, filter_invalid=True, ignore_UUIDs=False):
        """Returns an Nx2 matrix, 1st column is sequence, 2nd column is count"""

        if not self._has_UUIDs:
            ignore_UUIDs = True

        if not by_amino_acid:

            sequence_counts = {}
            
            if ignore_UUIDs:

                for sequence_UUID_count in self._sequence_UUID_counts:

                    sequence = sequence_UUID_count[0]

                    if sequence not in sequence_counts:
                        sequence_counts[sequence] = sequence_UUID_count[2]
                    else:
                        sequence_counts[sequence] += sequence_UUID_count[2]
                        
            else:
                for sequence_UUID_count in self._sequence_UUID_counts:

                    sequence = sequence_UUID_count[0]

                    if sequence not in sequence_counts:
                        sequence_counts[sequence] = 1
                    else:
                        sequence_counts[sequence] += 1

            filtered_sequence_counts = {}

            for sequence, sequence_count in sequence_counts.items():

                if sequence_count < count_threshold:
                    pass
                else:
                    filtered_sequence_counts[sequence] = sequence_count
                    
            return filtered_sequence_counts
        else:
            sequence_counts = {}
            
            if ignore_UUIDs:

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
            else:

                for sequence_UUID_count in self._sequence_UUID_counts:

                    sequence = sequence_UUID_count[0]

                    if (filter_invalid or by_amino_acid) and sequence.find('N') != -1:
                        continue

                    sequence_count = 1

                    if by_amino_acid:

                        sequence = DNA.translate_dna_single(sequence)

                        if filter_invalid and sequence.find("#") != -1:
                            continue

                    if sequence not in sequence_counts:
                        sequence_counts[sequence] = sequence_count
                    else:
                        sequence_counts[sequence] += sequence_count

        if count_threshold > 1:
            masked_sequence_counts = {}
            for sequence, sequence_count in sequence_counts.items():

                if sequence_count < count_threshold:
                    pass
                else:
                    masked_sequence_counts[sequence] = sequence_count

            return masked_sequence_counts
        else:
            return sequence_counts

    def eliminate_bias(self, predicted_bias_percentage=0.01, num_nucleotides_off=1, use_UIDs=False):

        if num_nucleotides_off != 1:
            raise NotImplementedError("Only support single nucleotide errors")

        if use_UIDs:
            raise NotImplementedError("Currently ignoring UIDs in eliminate bias")

        print("Getting sequence counts")
        # Get a dictionary of unique sequences and their respective counts
        sequence_count_dict = self.get_sequence_counts(by_amino_acid=False, count_threshold=0, filter_invalid=True)

        print("Converting sequence counts to sorted lists")
        sequence_counts = pandas.DataFrame.from_dict(sequence_count_dict, orient="index")
        sequence_counts.columns = ["count"]
        sequence_counts.sort_values(by="count", ascending=True, inplace=True)

        sequences = list(sequence_counts.index)
        counts = numpy.array(sequence_counts["count"].values)

        print("Building sequence trie")
        sequence_trie = Sequence_Trie(by_nucleotide=True)
        for sequence, count in zip(sequences, counts):
            sequence_trie.add(sequence, count)

        alphabet = set(DNA.get_nucleotides())

        for sequence_index, sequence in enumerate(sequences):

            if sequence_index % 10000 == 0:
                print("Analyzing sequence %i" % sequence_index)

            count = counts[sequence_index]

            min_parent_count = math.floor(count * (1 / predicted_bias_percentage))

            parent_counts = []

            # Find all possible parent sequences
            for index, character in enumerate(sequence):

                prefix = sequence[0:index]

                parent_node = sequence_trie.get_node(prefix)

                if parent_node is None:
                    continue

                for other_character in alphabet.difference(character):

                    postfix = other_character + sequence[index + 1:]

                    parent_count = parent_node.get_value(postfix)

                    if parent_count is None:
                        continue
                    elif parent_count < min_parent_count:
                        continue

                    parent_counts.append((prefix + postfix, parent_count))

            if len(parent_counts) > 0:

                parent_count_sum = 0

                for parent, parent_count in parent_counts:
                    parent_count_sum += parent_count

                for parent, parent_count in parent_counts:
                    count_to_add = parent_count / parent_count_sum * count
                    #             print("Adding '%.2f' to parent '%s' from '%s'" % (count_to_add, parent, sequence))
                    sequence_trie.add(parent, parent_count + count_to_add)
                sequence_trie.add(sequence, 0)

        self._sequence_UUID_counts = []

        for sequence in sequences:
            new_sequence_count = sequence_trie.get_value(sequence)
            if new_sequence_count > 0:
                self._sequence_UUID_counts.append((sequence, "", new_sequence_count))

    def sequence_comparable(sequence_to_compare,sequence_to_be_compared_to,threshold):
        """Returns true if two string are comparable given user defined threshold"""
        
        distance = 0
        
        for i in range(len(sequence_to_compare)):
            if sequence_to_compare[i] != sequence_to_be_compared_to[i]:
                distance += 1
            if distance > threshold:
                return False
        return True