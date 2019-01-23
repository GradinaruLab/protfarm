import random

from peseq.utils import DNA

from ..workspace import Workspace as ws
from ..workspace import Database as db
from . import Sequence_Library

def get_probability_of_unseen_sequence(library):

    alignment = ws.get_active_alignment()

    # Check the alignment statistics data to see if we've gotten calculated this before
    if "Unseen Sequence Probability" in alignment.statistics[library.id]:
        return alignment.statistics[library.id]["Unseen Sequence Probability"]

    if "Expected Number of Misreads" not in alignment.statistics[library.id]:
        raise Exception("Missing 'Expected Number of Misreads' from statistics. Align with an alignment method that generates ")

    sequence_library = Sequence_Library(library)

    num_expected_misreads = alignment.statistics[library.id]["Expected Number of Misreads"]
    sequence_counts = sequence_library.get_sequence_counts(by_amino_acid=False, count_threshold=0, filter_invalid=True)
    num_single_counts = 0

    for sequence, count in sequence_counts.items():
        if count == 1:
            num_single_counts += 1

    probability_of_misread_overlap = coverage.get_probability_of_single_misread_existing(library)
    print("num_expected_misreads: %.4f" % num_expected_misreads)
    print("num_single_counts: %i" % num_single_counts)
    print("probability_of_misread_overlap: %.4f" % probability_of_misread_overlap)
    unique_misreads = round(num_expected_misreads * (1-probability_of_misread_overlap))

    n_1 = num_single_counts - unique_misreads

    if n_1 <= num_single_counts * -10:
        raise Exception("Order of magnitude more expected unique misreads than there are single count sequences! This should be impossible")

    if n_1 <= 0:
        n_1 = 1

    probability_unseen = n_1 / alignment.statistics[library.id]["Number of Sequences"]

    alignment.set_statistic(library, "Unseen Sequence Probability", probability_unseen)

    return probability_unseen


def get_probability_of_single_misread_existing(library, num_samples = 10000):

    alignment = ws.get_active_alignment()

    template = db.get_template_by_id \
        (alignment.library_templates[library.id]).get_variant_template()

    sequence_library = Sequence_Library(library)

    sequence_counts = sequence_library.get_sequence_counts(by_amino_acid=False, count_threshold=0, filter_invalid=True)

    sequences = []

    for sequence, count in sequence_counts.items():
        for i in range(count):
            sequences.append(sequence)

    num_duplicates = 0
    sample_index = 0

    while sample_index < num_samples:
        sequence = random.sample(sequences, 1)[0]

        random_mutation_index = random.sample(range(len(template)), 1)[0]

        new_nucleotide = None

        while new_nucleotide == None or new_nucleotide == sequence[random_mutation_index]:
            new_nucleotide = random.sample(list(DNA.IUPAC_GRAMMAR_MAP[template[random_mutation_index]]), 1)[0]

        new_sequence = sequence[:random_mutation_index] + new_nucleotide + sequence[random_mutation_index +1:]

        sample_index += 1

        if new_sequence in sequence_counts:
            num_duplicates += 1

    return num_duplicates / sample_index


def get_coverage(analysis_set, by_amino_acid = True):

    num_included_sequences = 0

    variant_template = ''

    alignment = ws.get_active_alignment()

    # Get the templates for each analysis set
    for library_name in analysis_set.get_libraries():
        db_library = db.get_library(library_name)
        template_id = alignment.library_templates[db_library.id]
        template = db.get_template_by_id(template_id)
        if len(variant_template) == 0:
            variant_template = template.get_variant_template()
        elif variant_template != template.get_variant_template():
            raise Exception('Variant templates must match in an analysis set to do coverage analysis!')

    num_possible_sequences = 1

    if by_amino_acid:
        if len(variant_template) % 3 != 0:
            raise Exception('Can\'t analyze by amino acid when variant sequence isn\'t groups of 3!')
        variant_length = int(len(variant_template) / 3)

        # For now, assume all amino acids are possible. Should do something with IUPAC later
        for variant_index in range(0, variant_length):
            num_possible_sequences *= 20;
    else:
        for variant_index in range(0, len(variant_template)):
            num_possible_sequences *= len(DNA.IUPAC_GRAMMAR_MAP[variant_template[variant_index]])

    # We assume all sequences in the analysis set have been aligned against the template,
    # so they must match the template. So, all unique sequences are the included sequences

    included_sequences = set()

    for library_name, library in analysis_set.get_libraries().items():
        sequence_counts = library.get_sequence_counts(by_amino_acid, count_threshold = 0)

        for sequence, count in sequence_counts.items():
            included_sequences.add(sequence)

    return len(included_sequences), num_possible_sequences
