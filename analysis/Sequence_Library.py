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
