import math
import pandas
import numpy

from pepars.utils import DNA
from pepars.analysis import confidence
from pepars.analysis import Label_Type
from pepars.utils import Sequence_Trie

from . import coverage

from ..workspace import Workspace as ws
from ..workspace import Database as db

from .Sequence_Library import Sequence_Library


class Analysis_Set:

    def __init__(self):

        self.sequence_libraries = {}
        self.sequence_length = 0
        self._sample_sequence_counts = None
        self._sample_total_counts = None

    def add_library(self, library):

        if isinstance(library,str):
            library = db.get_library(library)

        if library.name in self.sequence_libraries:
            return

        sequence_library = Sequence_Library(library)
        self.sequence_length = sequence_library.get_sequence_length()
        self.sequence_libraries[library.name] = sequence_library

    # Preparing for renaming library -> sample
    def add_sample(self, sample_name):
        self.add_library(sample_name)

    def get_sequence_length(self, by_amino_acid=True):

        if by_amino_acid:
            return int(self.sequence_length / 3)
        else:
            return self.sequence_length

    def get_libraries(self):
        return self.sequence_libraries

    def load_sequence_counts(self, count_threshold=0):

        self._sample_sequence_counts = {}
        self._sample_total_counts = {}

        for sample_name in self.sequence_libraries:
            sample_counts_dict = self.get_libraries()[sample_name]\
                .get_sequence_counts(count_threshold=count_threshold)

            sequence_counts = pandas.DataFrame.from_dict(
                sample_counts_dict, orient="index")
            sequence_counts.columns = ["%s_count" % sample_name]
            sequence_counts.sort_index(inplace=True)
            self._sample_sequence_counts[sample_name] = sequence_counts
            self._sample_total_counts[sample_name] = sequence_counts.sum()

    # libraries_of_interest: list of identifiers for the library that you want the specificity of
    # libraries_to_compare: list of identifiers for the library that you want to compare the specificity against
    #
    # Returns: Nx2 matrix, 1st column is sequence, 2nd column is specificity score
    def get_specificity(self, library_of_interest_names,
        libraries_to_compare_names, by_amino_acid = True, count_threshold = 10,
        log_scale = False, zero_count_magic_number = None):

        alignment = ws.get_active_alignment()

        specificity_dict={}

        library_of_interest_total_count = 0

        if isinstance(library_of_interest_names, list):

            sequence_counts = {}

            library_of_interest_unseen_probabilities = []

            for library_of_interest_name in library_of_interest_names:
                library_of_interest = self.sequence_libraries[library_of_interest_name]
                library_of_interest_total_count += library_of_interest.get_total_count()
                library_of_interest_counts = library_of_interest.get_sequence_counts(by_amino_acid, count_threshold = 0)
                if zero_count_magic_number == None:
                    library_of_interest_unseen_probabilities.append(
                        coverage.get_probability_of_unseen_sequence(db.get_library(library_of_interest_name)))

                for sequence, count in library_of_interest_counts.items():
                    if sequence not in sequence_counts:
                        sequence_counts[sequence] = count
                    else:
                        sequence_counts[sequence] += count

            if zero_count_magic_number == None:
                library_of_interest_unseen_probability = sum(library_of_interest_unseen_probabilities) / len(library_of_interest_unseen_probabilities)
            else:
                library_of_interest_unseen_probability = zero_count_magic_number

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
            if zero_count_magic_number == None:
                library_of_interest_unseen_probability = coverage.get_probability_of_unseen_sequence(db.get_library(library_of_interest_names))
            else:
                library_of_interest_unseen_probability = zero_count_magic_number

            num_comparing_libraries = len(libraries_to_compare_names) + 1

        libraries_to_compare = []
        libraries_to_compare_total_counts = []
        libraries_to_compare_unseen_probabilities = []

        for library_to_compare_name in libraries_to_compare_names:
            library_to_compare = self.sequence_libraries[library_to_compare_name]
            libraries_to_compare_total_counts.append(library_to_compare.get_total_count())
            library_to_compare = library_to_compare.get_sequence_counts(by_amino_acid, count_threshold)
            if zero_count_magic_number == None:
                libraries_to_compare_unseen_probabilities.append(coverage.get_probability_of_unseen_sequence(db.get_library(library_to_compare_name)))
            else:
                libraries_to_compare_unseen_probabilities.append(zero_count_magic_number)
            libraries_to_compare.append(library_to_compare)

        num_sequences_that_exist_but_have_negative_specificity = 0
        num_sequences_that_dont_exist_in_libraries_to_compare = 0

        for key in sequence_counts:

            comparing_libraries_presence = 0

            exists_in_other_library = False

            for library_to_compare_index in range(len(libraries_to_compare)):
                library_to_compare = libraries_to_compare[library_to_compare_index]

                library_presence = 0
                if key in library_to_compare:
                    exists_in_other_library = True
                    #print key + " is in! Count: " + str(library_to_compare[key])
                    library_presence = library_to_compare[key] * 1.0 /libraries_to_compare_total_counts[library_to_compare_index]
                    comparing_libraries_presence += library_presence
                else:
                    library_presence = libraries_to_compare_unseen_probabilities[library_to_compare_index] /libraries_to_compare_total_counts[library_to_compare_index]

            library_of_interest_presence = 1.0*sequence_counts[key]/library_of_interest_total_count
            comparing_libraries_presence += library_of_interest_presence

            if not exists_in_other_library:
                num_sequences_that_dont_exist_in_libraries_to_compare += 1

            sequence_specificity = library_of_interest_presence/(comparing_libraries_presence/num_comparing_libraries)

            if log_scale:
                specificity_dict[key] = math.log10(sequence_specificity)
                if math.log10(sequence_specificity) < 0:
                    num_sequences_that_exist_but_have_negative_specificity += 1
            else:
                specificity_dict[key] = sequence_specificity

        print('Of the %i sequences that exist, %i of them have negative specificity, (%f%%)' % (len(sequence_counts), num_sequences_that_exist_but_have_negative_specificity, num_sequences_that_exist_but_have_negative_specificity/len(sequence_counts)))
        print('And %i of them don\'t exist in the libraries to compare to' % num_sequences_that_dont_exist_in_libraries_to_compare)

        for library_to_compare_index in range(len(libraries_to_compare)):
            library_to_compare = libraries_to_compare[library_to_compare_index]
            for key in library_to_compare:
                if key not in specificity_dict:

                    library_of_interest_presence = library_of_interest_unseen_probability/library_of_interest_total_count
                    comparing_libraries_presence = library_to_compare[key] / libraries_to_compare_total_counts[library_to_compare_index]

                    for other_library_to_compare_index in range(library_to_compare_index + 1,len(libraries_to_compare)):
                        other_library_to_compare = libraries_to_compare[other_library_to_compare_index]
                        if key in other_library_to_compare:
                            comparing_libraries_presence += other_library_to_compare[key] / libraries_to_compare_total_counts[other_library_to_compare_index]
                        else:
                            comparing_libraries_presence += libraries_to_compare_unseen_probabilities[library_to_compare_index] / libraries_to_compare_total_counts[other_library_to_compare_index]

                    comparing_libraries_presence += library_of_interest_presence
                    sequence_specificity = library_of_interest_presence/(comparing_libraries_presence/num_comparing_libraries)

                    if log_scale:
                        specificity_dict[key] = math.log10(sequence_specificity)
                    else:
                        specificity_dict[key] = sequence_specificity

        return specificity_dict

    def get_enrichment(self, library_of_interest_names, starting_library_names,
        by_amino_acid = True, count_threshold = 10, Log_Scale=True,
        zero_count_magic_number = None, include_zero_count=True,
        filter_invalid = True, count_threshold_starting_library = None,
        count_threshold_library_of_interest = None,
        include_zero_count_starting_library = None,
        include_zero_count_library_of_interest = None,
        zero_count_magic_number_starting_library = None,
        zero_count_magic_number_library_of_interest = None):

        if zero_count_magic_number_starting_library == None:
            zero_count_magic_number_starting_library = zero_count_magic_number
        if include_zero_count_starting_library == None:
            include_zero_count_starting_library = include_zero_count

        if zero_count_magic_number_library_of_interest == None:
            zero_count_magic_number_library_of_interest = zero_count_magic_number
        if include_zero_count_library_of_interest == None:
            include_zero_count_library_of_interest = include_zero_count

        if isinstance(library_of_interest_names, list):

            library_of_interest_total_count = 0
            sequence_counts = {}
            library_of_interest_unseen_probabilities = []

            for library_of_interest_name in library_of_interest_names:
                library_of_interest = self.sequence_libraries[library_of_interest_name]
                library_of_interest_total_count += library_of_interest.get_total_count()
                library_of_interest_counts = library_of_interest.get_sequence_counts(by_amino_acid=by_amino_acid, count_threshold = 0, filter_invalid = filter_invalid)

                if zero_count_magic_number_library_of_interest == None and include_zero_count_library_of_interest:
                    library_of_interest_unseen_probabilities.append(coverage.get_probability_of_unseen_sequence(db.get_library(library_of_interest_name)))

                for sequence, count in library_of_interest_counts.items():
                    if sequence not in sequence_counts:
                        sequence_counts[sequence] = count
                    else:
                        sequence_counts[sequence] += count

            if zero_count_magic_number_library_of_interest == None and include_zero_count_library_of_interest:
                library_of_interest_unseen_probability = sum(library_of_interest_unseen_probabilities) / len(library_of_interest_unseen_probabilities)
            else:
                library_of_interest_unseen_probability = zero_count_magic_number_library_of_interest

            library_of_interest = sequence_counts
        else:
            library_of_interest = self.sequence_libraries[library_of_interest_names]
            print("Getting total count of %s" % library_of_interest_names)
            library_of_interest_total_count = library_of_interest.get_total_count()
            print("Got total count of %s" % library_of_interest_names)

            print("Getting probability of unseen sequence for %s" % library_of_interest_names)
            if zero_count_magic_number_library_of_interest == None and include_zero_count_library_of_interest:
                library_of_interest_unseen_probability = coverage.get_probability_of_unseen_sequence(db.get_library(library_of_interest_names))
            else:
                library_of_interest_unseen_probability = zero_count_magic_number_library_of_interest
            print("Got probability of unseen sequence for %s" % library_of_interest_names)

            print("Getting sequence counts of %s" % library_of_interest_names)
            library_of_interest = library_of_interest.get_sequence_counts(by_amino_acid=by_amino_acid, count_threshold = 0, filter_invalid = filter_invalid)
            print("Got sequence counts of %s" % library_of_interest_names)

        if isinstance(starting_library_names, list):

            aggregate_starting_library = {}
            starting_library_unseen_probabilities = []

            for starting_library_name in starting_library_names:
                starting_library = self.sequence_libraries[starting_library_name]
                starting_library_total_count = starting_library.get_total_count()
                starting_library = starting_library.get_sequence_counts(by_amino_aci=by_amino_acid, count_threshold = 0, filter_invalid = filter_invalid)

                if zero_count_magic_number_starting_library == None and include_zero_count_starting_library:
                    starting_library_unseen_probabilities.append(coverage.get_probability_of_unseen_sequence(db.get_library(starting_library_name)))

                for sequence, sequence_count in starting_library.items():

                    if sequence not in aggregate_starting_library:
                        aggregate_starting_library[sequence] = sequence_count
                    else:
                        aggregate_starting_library[sequence] += sequence_count

            if zero_count_magic_number_starting_library == None and include_zero_count_starting_library:
                starting_library_unseen_probability = sum(starting_library_unseen_probabilities) / len(starting_library_unseen_probabilities)
            else:
                starting_library_unseen_probability = zero_count_magic_number_starting_library

            starting_library = aggregate_starting_library

        else:
            starting_library = self.sequence_libraries[starting_library_names]
            starting_library_total_count = starting_library.get_total_count()

            if zero_count_magic_number_starting_library == None and include_zero_count_starting_library:
                starting_library_unseen_probability = coverage.get_probability_of_unseen_sequence(db.get_library(starting_library_names))
            else:
                starting_library_unseen_probability = zero_count_magic_number_starting_library

            starting_library = starting_library.get_sequence_counts(by_amino_acid=by_amino_acid, count_threshold = 0, filter_invalid = filter_invalid)

        # If the user has specified separate thresholds for the starting
        # library and library of interest, we process them separately
        if count_threshold_starting_library != None or count_threshold_library_of_interest != None:
            if count_threshold_starting_library == None:
                count_threshold_starting_library = count_threshold
            if count_threshold_library_of_interest == None:
                count_threshold_library_of_interest = count_threshold

            below_threshold_sequences = set()

            for sequence, count in starting_library.items():
                if count < count_threshold_starting_library:
                    below_threshold_sequences.add(sequence)

            for sequence in below_threshold_sequences:
                del(starting_library[sequence])

            below_threshold_sequences = set()

            for sequence, count in library_of_interest.items():
                if count < count_threshold_library_of_interest:
                    below_threshold_sequences.add(sequence)

            for sequence in below_threshold_sequences:
                del(library_of_interest[sequence])

        else:
            below_threshold_sequences = set()

            # Find all sequences in the library of interest that are below threshold
            for sequence, count in library_of_interest.items():
                if sequence in starting_library:
                    count += starting_library[sequence]

                if count < count_threshold:
                    below_threshold_sequences.add(sequence)

            # Find all sequences in the starting library that are below threshold
            for sequence, count in starting_library.items():
                if sequence in library_of_interest:
                    count += library_of_interest[sequence]

                if count < count_threshold:
                    below_threshold_sequences.add(sequence)

            # Remove below threshold sequences from both libraries
            for sequence in below_threshold_sequences:
                if sequence in library_of_interest:
                    del(library_of_interest[sequence])
                if sequence in starting_library:
                    del(starting_library[sequence])

        for sequence in starting_library:
            if sequence not in library_of_interest and include_zero_count_library_of_interest:
                library_of_interest_total_count += library_of_interest_unseen_probability

        for sequence in library_of_interest:
            if sequence not in starting_library and include_zero_count_starting_library:
                starting_library_total_count += starting_library_unseen_probability

        enrichment_dict = {}
        for sequence in starting_library:

            if sequence not in library_of_interest:
                if include_zero_count_library_of_interest:
                    library_of_interest[sequence] = library_of_interest_unseen_probability
                else:
                    continue

            if (Log_Scale):
                fold_enrichment = (library_of_interest[sequence]* 1.0 / library_of_interest_total_count) / (starting_library[sequence] * 1.0/ starting_library_total_count)
                enrichment_dict[sequence] = math.log10(fold_enrichment)
            else:
                fold_enrichment = (library_of_interest[sequence]* 1.0 / library_of_interest_total_count) / (starting_library[sequence] * 1.0/ starting_library_total_count)
                enrichment_dict[sequence] = fold_enrichment

        # For the sequences that don't exist in the starting library
        if (include_zero_count_starting_library):
            for sequence in library_of_interest:

                if sequence in starting_library:
                    continue

                fold_enrichment = (library_of_interest[sequence] * 1.0 / library_of_interest_total_count) / (starting_library_unseen_probability / starting_library_total_count)

                if Log_Scale:
                    enrichment_dict[sequence] = math.log10(fold_enrichment)
                else:
                    enrichment_dict[sequence] = fold_enrichment

        return enrichment_dict

    def get_sequence_weights(self, library_of_interest_names, starting_library_names,
        zero_count_magic_number = None,
        zero_count_magic_number_starting_library = None,
        zero_count_magic_number_library_of_interest = None,
        filter_invalid = True,
        by_amino_acid = True):

        if zero_count_magic_number_library_of_interest is None:
            zero_count_magic_number_library_of_interest = zero_count_magic_number

        if zero_count_magic_number_starting_library is None:
            zero_count_magic_number_starting_library = zero_count_magic_number

        if isinstance(library_of_interest_names, list):

            libraries_of_interest_sequence_counts = {}
            library_of_interest_unseen_probabilities = []

            for library_of_interest_name in library_of_interest_names:
                library_of_interest = self.sequence_libraries[library_of_interest_name]
                library_of_interest_sequence_counts = library_of_interest.get_sequence_counts(by_amino_acid=by_amino_acid, count_threshold = 0, filter_invalid = filter_invalid)

                if zero_count_magic_number_library_of_interest == None:
                    library_of_interest_unseen_probabilities.append(coverage.get_probability_of_unseen_sequence(db.get_library(library_of_interest_name)))

                for sequence, count in library_of_interest_sequence_counts.items():
                    if sequence not in libraries_of_interest_sequence_counts:
                        libraries_of_interest_sequence_counts[sequence] = count
                    else:
                        libraries_of_interest_sequence_counts[sequence] += count

            if zero_count_magic_number_library_of_interest == None:
                library_of_interest_unseen_probability = sum(library_of_interest_unseen_probabilities) / len(library_of_interest_unseen_probabilities)
            else:
                library_of_interest_unseen_probability = zero_count_magic_number_library_of_interest

        else:
            library_of_interest = self.sequence_libraries[library_of_interest_names]

            if zero_count_magic_number_library_of_interest == None:
                library_of_interest_unseen_probability = coverage.get_probability_of_unseen_sequence(db.get_library(library_of_interest_names))
            else:
                library_of_interest_unseen_probability = zero_count_magic_number_library_of_interest

            libraries_of_interest_sequence_counts = library_of_interest.get_sequence_counts(by_amino_acid=by_amino_acid, count_threshold = 0, filter_invalid = filter_invalid)

        if isinstance(starting_library_names, list):

            starting_libraries_sequence_counts = {}
            starting_library_unseen_probabilities = []

            for starting_library_name in starting_library_names:
                starting_library = self.sequence_libraries[starting_library_name]
                starting_library_sequence_counts = starting_library.get_sequence_counts(by_amino_acid=by_amino_acid, count_threshold = 0, filter_invalid = filter_invalid)

                if zero_count_magic_number_starting_library == None:
                    starting_library_unseen_probabilities.append(coverage.get_probability_of_unseen_sequence(db.get_library(starting_library_name)))

                for sequence, count in starting_library_sequence_counts.items():
                    if sequence not in starting_libraries_sequence_counts:
                        starting_libraries_sequence_counts[sequence] = count
                    else:
                        starting_libraries_sequence_counts[sequence] += count

            if zero_count_magic_number_starting_library == None:
                starting_library_unseen_probability = sum(starting_library_unseen_probabilities) / len(starting_library_unseen_probabilities)
            else:
                starting_library_unseen_probability = zero_count_magic_number_starting_library
        else:
            starting_library = self.sequence_libraries[starting_library_names]

            if zero_count_magic_number_starting_library == None:
                starting_library_unseen_probability = coverage.get_probability_of_unseen_sequence(db.get_library(starting_library_names))
            else:
                starting_library_unseen_probability = zero_count_magic_number_starting_library

            starting_libraries_sequence_counts = starting_library.get_sequence_counts(by_amino_acid=by_amino_acid, count_threshold = 0, filter_invalid = filter_invalid)

        sequence_weights = {}

        for sequence, count in libraries_of_interest_sequence_counts.items():
            sequence_weights[sequence] = count

            if sequence not in starting_libraries_sequence_counts and starting_library_unseen_probability is not None:
                sequence_weights[sequence] += starting_library_unseen_probability

        for sequence, count in starting_libraries_sequence_counts.items():

            if sequence not in sequence_weights:
                sequence_weights[sequence] = count
            else:
                sequence_weights[sequence] += count

            if sequence not in libraries_of_interest_sequence_counts and library_of_interest_unseen_probability is not None:
                sequence_weights[sequence] += library_of_interest_unseen_probability

        for sequence, count in sequence_weights.items():
            sequence_weights[sequence] = math.log2(count)

        return sequence_weights

    def export_enrichment(self, filename, starting_library_name, \
        by_amino_acid = False, count_threshold = 0, log_scale = False, filter_invalid = True,
        include_zero_count = True, zero_count_magic_number = 0.9):

        num_sequence_libraries = len(self.sequence_libraries)

        cumulative_counts = {}
        cumulative_enrichments = {}

        library_index = 0
        starting_library_index = None

        header_row = ['Sequence']
        if not by_amino_acid:
            header_row.append('Amino Acid')
        library_names = []

        for library_name, library in self.sequence_libraries.items():

            header_row.append(library_name)
            library_names.append(library_name)

            if library_name != starting_library_name:
                header_row.append(library_name + ' enrichment')

                if isinstance(starting_library_name, dict):
                    fold_enrichments = self.get_enrichment(library_name, starting_library_name[library_name], by_amino_acid, count_threshold = 0, Log_Scale = log_scale, include_zero_count = include_zero_count, zero_count_magic_number = zero_count_magic_number, filter_invalid = filter_invalid)
                else:
                    fold_enrichments = self.get_enrichment(library_name, starting_library_name, by_amino_acid, count_threshold = 0, Log_Scale = log_scale, include_zero_count = include_zero_count,  zero_count_magic_number = zero_count_magic_number, filter_invalid = filter_invalid)
            else:
                fold_enrichments = {}

            library_counts = library.get_sequence_counts(by_amino_acid, count_threshold = 0, filter_invalid = filter_invalid)

            # Look through all the counts we got
            for sequence, sequence_count in library_counts.items():

                if sequence not in cumulative_counts:
                    cumulative_counts[sequence] = {}

                cumulative_counts[sequence][library_name] = sequence_count

            # And look through all the enrichments we got
            for sequence, enrichment in fold_enrichments.items():

                if sequence not in cumulative_enrichments:
                    cumulative_enrichments[sequence] = {}

                cumulative_enrichments[sequence][library_name] = enrichment

        data = []

        num_sequences_filtered = 0

        for sequence, sequence_counts in cumulative_counts.items():

            sequence_row = [sequence]

            if not by_amino_acid:
                sequence_row.append(DNA.translate_DNA_to_AA(sequence))

            if sum(sequence_counts.values()) < count_threshold:
                num_sequences_filtered += 1
                continue

            library_index = 0

            for library_name in library_names:

                if library_name not in sequence_counts:
                    sequence_row.append(0)
                else:
                    sequence_row.append(sequence_counts[library_name])

                if library_name != starting_library_name:
                    sequence_row.append(cumulative_enrichments[sequence][library_name])

            data.append(sequence_row)

        ws.export_csv(filename, header_row, data)


    def export_enrichment_specificity(self, filename, starting_library_name, \
        libraries_to_compare_names, by_amino_acid = False, \
        count_threshold = 0, log_scale = True, include_zero_count = True, zero_count_magic_number = None):

        if not starting_library_name:
            calculate_enrichment = False
        else:
            calculate_enrichment = True

        if not libraries_to_compare_names:
            calculate_specificity = False
        else:
            calculate_specificity = True

        num_sequence_libraries = len(self.sequence_libraries)

        cumulative_counts = {}
        cumulative_enrichments = {}
        cumulative_specificities = {}

        library_index = 0

        header_row = ['Sequence']
        if not by_amino_acid:
            header_row.append('Amino Acid')

        for library_name, library in self.sequence_libraries.items():

            print("Getting enrichment for %s" % library_name)

            if library_name == starting_library_name:
                library_count_threshold = count_threshold
            else:
                library_count_threshold = 0

            header_row.append(library_name + ' count')
            if calculate_enrichment:
                header_row.append(library_name + ' enrichment')
                fold_enrichments = self.get_enrichment(library_name,
                                                       starting_library_name,
                                                       by_amino_acid=by_amino_acid,
                                                       count_threshold=library_count_threshold,
                                                       Log_Scale=log_scale,
                                                       include_zero_count=include_zero_count,
                                                       zero_count_magic_number=zero_count_magic_number)

            print("Got enrichment for %s" % library_name)

            if calculate_specificity:
                header_row.append(library_name + ' specificity')
                specificities = self.get_specificity(library_name, libraries_to_compare_names, by_amino_acid = by_amino_acid, count_threshold = 0, log_scale = log_scale, zero_count_magic_number = zero_count_magic_number)

            print("Getting library counts for %s" % library_name)

            library_counts = library.get_sequence_counts(by_amino_acid, count_threshold = 0, filter_invalid = False)

            print("Got library counts for %s" % library_name)

            for sequence, sequence_count in library_counts.items():

                if sequence not in cumulative_counts:
                    cumulative_counts[sequence] = [0] * num_sequence_libraries
                    cumulative_enrichments[sequence] = [0] * num_sequence_libraries
                    cumulative_specificities[sequence] = [0] * num_sequence_libraries

                if calculate_enrichment and sequence in fold_enrichments:
                    cumulative_enrichments[sequence][library_index] = fold_enrichments[sequence]

                if calculate_specificity and sequence in specificities:
                    cumulative_specificities[sequence][library_index] = specificities[sequence]

                cumulative_counts[sequence][library_index] = sequence_count

            library_index += 1

        data = []

        for sequence, sequence_counts in sorted(cumulative_counts.items(), key=lambda kv: sum(kv[1]), reverse=True):

            sequence_row = [sequence]

            if not by_amino_acid:
                sequence_row.append(DNA.translate_DNA_to_AA(sequence))

            if sum(sequence_counts) < count_threshold:
                continue

            library_index = 0

            if calculate_enrichment:
                fold_enrichments = cumulative_enrichments[sequence]
            if calculate_specificity:
                specificities = cumulative_specificities[sequence]

            for sequence_count in sequence_counts:
                sequence_row.append(sequence_count)
                if calculate_enrichment:
                    sequence_row.append(fold_enrichments[library_index])
                if calculate_specificity:
                    sequence_row.append(specificities[library_index])
                library_index += 1

            data.append(sequence_row)

        ws.export_csv(filename, header_row, data)

    def get_sequence_labels_confidences(
            self,
            target_samples,
            starting_sample=None,
            off_target_samples=None,
            confidence_metric=confidence.Confidence_Metric.LOG_PROBABILITY,
            label_type=Label_Type.ENRICHMENT,
            renormalize_totals=False,
            target_sample_weights=None,
            off_target_sample_weights=None,
            filter_sequences=None,
            standardize=False,
            exclude_sequences=None):

        if self._sample_sequence_counts is None:
            raise EnvironmentError("Sample counts not loaded, must call "
                                   "load_sample_counts before calling "
                                   "get_sequence_labels_confidences")

        if isinstance(target_samples, str):
            target_samples = [target_samples]

        if off_target_samples is not None and len(off_target_samples) == 0:
            off_target_samples = None

        if isinstance(off_target_samples, str):
            off_target_samples = [off_target_samples]

        if starting_sample is None and label_type in \
                [Label_Type.ENRICHMENT, Label_Type.ENRICHMENT_SPECIFICITY]:
            raise ValueError("Must have starting sample unless predicting on "
                             "counts")

        if off_target_samples is None and label_type in \
                [Label_Type.SPECIFICITY, Label_Type.ENRICHMENT_SPECIFICITY]:
            raise ValueError("Must have off target samples if looking at "
                             "specificity")

        sample_names = set()

        for sample in target_samples:
            sample_names.add(sample)

        if off_target_samples is not None:
            for sample in off_target_samples:
                sample_names.add(sample)

        if starting_sample is not None:
            sample_names.add(starting_sample)

        count_data_frames = []

        for sample_name in sample_names:
            count_data_frames.append(self._sample_sequence_counts[sample_name])

        sequence_counts = pandas.concat(count_data_frames, axis=1, sort=True)

        if exclude_sequences is not None:
            for sequence in exclude_sequences:
                if sequence not in sequence_counts.index:
                    continue
                sequence_counts.drop(sequence, inplace=True)

        for sample_name in sample_names:
            sequence_counts["%s_count" % sample_name] = \
                sequence_counts["%s_count" % sample_name].fillna(
                    coverage.get_probability_of_unseen_sequence(
                        db.get_library(sample_name)
                    )
                )

        sequence_counts["Confidence"] = confidence.get_sequence_confidences(
            sequence_counts.values,
            confidence_metric=confidence_metric
        )

        labels_confidences = pandas.DataFrame(sequence_counts["Confidence"])

        if label_type in [Label_Type.LOG_COUNTS, Label_Type.RAW]:
            if off_target_samples is not None:
                raise ValueError("Should not specify off target samples when"
                                 "looking at counts")
            count_column_names = ["%s_count" % name for name in sample_names]
            labels_confidences["Label"] = sequence_counts[count_column_names].sum(axis=1)

            if label_type == Label_Type.LOG_COUNTS:
                labels_confidences["Label"] = numpy.log(labels_confidences["Label"])

        for sample_name in sample_names:
            if renormalize_totals:
                sequence_counts["%s_probability" % sample_name] = \
                    sequence_counts["%s_count" % sample_name] / \
                    sequence_counts["%s_count" % sample_name].sum()
            else:
                sequence_counts["%s_probability" % sample_name] = \
                    sequence_counts["%s_count" % sample_name] / \
                    self._sample_total_counts[sample_name].values

        target_count_column_names = \
            ["%s_probability" % name for name in target_samples]

        sequence_counts["target_probability"] = \
            sequence_counts[target_count_column_names].mean(axis=1)

        if label_type in \
                [Label_Type.ENRICHMENT, Label_Type.ENRICHMENT_SPECIFICITY]:
            sequence_counts["enrichment"] = \
                numpy.log(sequence_counts["target_probability"] /
                          sequence_counts["%s_probability" % starting_sample])

        if label_type in \
                [Label_Type.SPECIFICITY, Label_Type.ENRICHMENT_SPECIFICITY]:

            off_target_count_column_names = \
                ["%s_probability" % name for name in off_target_samples]

            sequence_counts["off_target_probability"] = \
                sequence_counts[off_target_count_column_names].mean(axis=1)

            sequence_counts["specificity"] = \
                numpy.log(sequence_counts["target_probability"] /
                          sequence_counts["off_target_probability"])

        if label_type == Label_Type.ENRICHMENT:
            labels_confidences["Label"] = sequence_counts["enrichment"]
        elif label_type == Label_Type.SPECIFICITY:
            labels_confidences["Label"] = sequence_counts["specificity"]
        elif label_type == Label_Type.ENRICHMENT_SPECIFICITY:
            sequence_counts["enrichment"] = sequence_counts["enrichment"] - \
                sequence_counts["enrichment"].min()
            sequence_counts["enrichment"] = sequence_counts["enrichment"] / \
                sequence_counts["enrichment"].max()
            sequence_counts["specificity"] = sequence_counts["specificity"] - \
                sequence_counts["specificity"].min()
            sequence_counts["specificity"] = sequence_counts["specificity"] / \
                sequence_counts["specificity"].max()
            labels_confidences["Label"] = \
                sequence_counts["enrichment"] * sequence_counts["specificity"]

        if standardize:
            labels_confidences["Label"] -= \
                labels_confidences["Label"].mean()
            labels_confidences["Label"] /= \
                labels_confidences["Label"].std()
        else:

            labels_confidences["Label"] = labels_confidences["Label"] - \
                labels_confidences["Label"].min()
            labels_confidences["Label"] = labels_confidences["Label"] / \
                labels_confidences["Label"].max()

        for sample_name in sample_names:
            labels_confidences["%s_count" % sample_name] = \
                sequence_counts["%s_count" % sample_name]

        return labels_confidences
