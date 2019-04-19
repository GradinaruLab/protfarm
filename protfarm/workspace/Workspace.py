import os
from . import Library
from . import Database as db
from . import FASTQ_File
from . import Alignment

from pepars.fileio import csv_wrapper
from pepars.utils import FASTQ_File
from pepars.utils import FASTQ_File_Set
from pepars.alignment import Aligner
from pepars.alignment.Perfect_Match_Aligner import Perfect_Match_Aligner
from pepars.alignment.Bowtie_Aligner import Bowtie_Aligner


def get_fastq_file_names():

    raw_data_directory = get_full_path(raw_data_subdirectory)
    mkdir_if_not_exists(raw_data_directory)

    fastq_files = [file for file in os.listdir(raw_data_directory) \
        if file.endswith('.fastq') or file.endswith('fastq.gz')]

    # Truncate the .gz since the user doesn't need to care whether it's a .gz or not
    for i in range(0, len(fastq_files)):
        if fastq_files[i].endswith('.gz'):
            fastq_files[i] = fastq_files[i][:-3]

    fastq_files.sort()

    return fastq_files


def write_sequence_file(library, alignment, sequence_uuid_counts):

    file_name = get_alignment_file_name(library, alignment)

    csv_wrapper.write_csv_file(file_name, ['Sequence', 'UUID', 'Count'], \
        sequence_uuid_counts)


def get_alignment_file_name(library, alignment):

    aligned_directory = get_full_path(aligned_subdirectory)
    mkdir_if_not_exists(aligned_directory)

    file_name = aligned_directory + "/" + \
        str(library.id) + "_" + str(alignment.id) + ".csv"

    return file_name


def alignment_exists(library, alignment):

    alignment_file_name = get_alignment_file_name(library, alignment)

    if os.path.isfile(alignment_file_name) and \
        library.id in alignment.statistics:
        return True

    return False


def remove_library_alignments(library):

    alignments = db.get_alignments()

    for alignment in alignments:
        file_name = get_alignment_file_name(library, alignment)
        alignment.remove_library(library)
    try:
        os.remove(file_name)
    except OSError:
        pass


def cleanup(workspace_path):
    """Clean up the workspace of any temporary files from an unclean shutdown"""

    raw_data_directory = get_full_path(raw_data_subdirectory)
    mkdir_if_not_exists(raw_data_directory)

    fastq_files = [file for file in os.listdir(raw_data_directory) \
        if file.endswith('.fastq')]

    compressed_fastq_files = [file for file in os.listdir(raw_data_directory) \
        if file.endswith('fastq.gz')]

    # Check if there are any uncompressed files that match a compressed file. If
    # so we delete the raw data. Otherwise we add it to our list

    for compressed_fastq_file in compressed_fastq_files:

        fastq_file_name = compressed_fastq_file[:-3]

        already_exists = False

        for fastq_file_index in range(0, len(fastq_files)):
            if fastq_files[fastq_file_index] == fastq_file_name:
                fastq_file_path = os.path.join(raw_data_directory, \
                    fastq_file_name)
                os.remove(fastq_file_path)
                already_exists = True
                break

        if not already_exists:
            fastq_files.append(fastq_file_name)

    fastq_files.sort()

    return fastq_files


def set_data_path(new_data_path):

    global data_path

    data_path = new_data_path


def set_experiment(new_experiment_name):

    global data_path

    if data_path is None and "VIRUS_FARM_LIBRARY_PATH" not in os.environ:
        raise EnvironmentError("VIRUS_FARM_LIBRARY_PATH not defined as an "
                               "environment variable.")

    
    if data_path is None:
        data_path = os.environ["VIRUS_FARM_LIBRARY_PATH"]

    global experiment_name

    experiment_name = new_experiment_name

    new_workspace_path = os.path.join(data_path, experiment_name)

    set_workspace_path(new_workspace_path)


def seed_libraries(library_illumina_project_map):

    FASTQ_file_names = get_fastq_file_names()
    illumina_project_files = {}

    library_names = library_illumina_project_map.keys()

    for FASTQ_file_name in FASTQ_file_names:
        parts = FASTQ_file_name.split("_")
        illumina_project_name = parts[0]
        if illumina_project_name not in illumina_project_files:
            illumina_project_files[illumina_project_name] = []

        illumina_project_files[illumina_project_name].append(FASTQ_file_name)

    for library_name in library_names:

        library = Library(library_name)

        illumina_project_name = library_illumina_project_map[library_name]

        if illumina_project_name not in illumina_project_files:
            raise ValueError("Missing '%s' from FASTQ files"
                             % illumina_project_name)

        for FASTQ_file_name in illumina_project_files[illumina_project_name]:
            library.add_file(FASTQ_file_name)


def set_workspace_path(new_workspace_path):
    """Set the current workspace path. This is where the db and all other
    files are located"""

    global workspace_path

    workspace_path = new_workspace_path

    try:
        mkdir_if_not_exists(os.path.join(workspace_path, raw_data_subdirectory))
    except:
        pass

    try:
        mkdir_if_not_exists(os.path.join(workspace_path, aligned_subdirectory))
    except:
        pass

    try:
        mkdir_if_not_exists(os.path.join(workspace_path, export_subdirectory))
    except:
        pass

    cleanup(workspace_path)

    db.load_database(workspace_path)

    global active_alignment

    try:
        active_alignment = db.get_active_alignment()
    except:
        active_alignment = None

    current_FASTQ_file_names = get_fastq_file_names()

    database_FASTQ_files = db.get_FASTQ_files()

    for FASTQ_file in database_FASTQ_files:
        if FASTQ_file.name not in current_FASTQ_file_names:
            continue#raise Exception("Missing previously existing '%s' FASTQ file!" % FASTQ_file.name)
        current_FASTQ_file_names.remove(FASTQ_file.name)

    for FASTQ_file_name in current_FASTQ_file_names:
        new_FASTQ_file = FASTQ_File(FASTQ_file_name)


def get_raw_data_path(child_path):
    return workspace_path + "/" + raw_data_subdirectory + "/" + child_path


def get_fastq_file(fastq_file_name):

    fastq_file_path = get_raw_data_path(fastq_file_name)

    if not os.path.exists(fastq_file_path):
        fastq_file_path = fastq_file_path + ".gz"

    if not os.path.exists(fastq_file_path):
        raise EnvironmentError("FASTQ file '%s' doesn't exist" %
                               fastq_file_name)

    fastq_file = FASTQ_File(fastq_file_path)

    fastq_files[fastq_file_name] = fastq_file

    return fastq_file


def get_full_path(child_path):
    return os.path.join(workspace_path, child_path)


def mkdir_if_not_exists(dir):

    if not os.path.isdir(dir):
        os.makedirs(dir)


def align_all(callback):

    global alignment_progress_string
    global alignment_progress_callback

    alignment_progress_callback = callback

    alignments = db.get_alignments()

    num_alignments = 0

    for alignment in alignments:
        if alignment.method not in [Perfect_Match_Aligner.__name__,
                Bowtie_Aligner.__name__]:
            raise Exception('Invalid alignment method, \'' + alignment.method \
                + '\', detected.')

        for library_id, template_id in alignment.library_templates.items():

            library = db.get_library_by_id(library_id)

            if not alignment_exists(library, alignment):
                num_alignments += 1

    alignment_number = 1

    for alignment in alignments:

        validate_alignment(alignment)

        for library_id, template_id in alignment.library_templates.items():

            library = db.get_library_by_id(library_id)
            template = db.get_template_by_id(template_id)

            if template.reverse_complement_template_id is not None:
                reverse_complement_template = \
                    db.get_template_by_id(
                        template.reverse_complement_template_id)
            else:
                reverse_complement_template = None

            if alignment_exists(library, alignment):
                continue

            alignment_progress_string = 'Performing alignment ' + \
                str(alignment_number) + '/' + str(num_alignments) + ': \'' \
                + library.name + '\' with \'' + alignment.method + '\''

            callback(alignment_progress_string)

            alignment_number += 1

            if alignment.method == Perfect_Match_Aligner.__name__:
                aligner = Perfect_Match_Aligner()
            elif alignment.method == Bowtie_Aligner.__name__:
                aligner = Bowtie_Aligner()
            else:
                continue

            read_1_FASTQ_files = []
            FASTQ_file_pairs = []
            FASTQ_file_sets = []

            for FASTQ_file_name in library.fastq_files:

                # Check to see if this is a read 1 file
                if FASTQ_file_name.find("_R1_") != -1:
                    read_1_FASTQ_files.append(FASTQ_file_name)

            for FASTQ_file_name in library.fastq_files:

                # Check to see if this is a read 2 file
                if FASTQ_file_name.find("_R2_") != -1:
                    read_1_FASTQ_files.append(FASTQ_file_name)
                    read_1_equivalent_name = \
                        FASTQ_file_name.replace("_R2_", "_R1_")

                    if read_1_equivalent_name in read_1_FASTQ_files:
                        FASTQ_file_pairs.append((read_1_equivalent_name,
                                                 FASTQ_file_name))

            if len(FASTQ_file_pairs) > 0:
                if len(FASTQ_file_pairs) != len(read_1_FASTQ_files):
                    raise ValueError("Must have equal numbers of read 1 and "
                                     "read 2 files")

                for FASTQ_file_pair in FASTQ_file_pairs:
                    read_1_file_path = get_raw_data_path(FASTQ_file_pair[0])
                    read_2_file_path = get_raw_data_path(FASTQ_file_pair[1])
                    FASTQ_file_sets.append(FASTQ_File_Set([read_1_file_path,
                                                           read_2_file_path]))
            else:
                for FASTQ_file_name in read_1_FASTQ_files:
                    read_1_file_path = get_raw_data_path(FASTQ_file_name)
                    FASTQ_file_sets.append(FASTQ_File_Set([read_1_file_path]))

            sequence_uuid_counts_array, statistics = \
                aligner.align(
                    template,
                    FASTQ_file_sets,
                    alignment.parameters,
                    reverse_complement_template=reverse_complement_template,
                    progress_callback=update_library_alignment_progress)

            write_sequence_file(library, alignment, sequence_uuid_counts_array)

            Alignment.add_statistics(library, statistics)


def validate_alignment(alignment):

    previous_num_variants = -1

    for library, template_id in alignment.library_templates.items():

        template = db.get_template_by_id(template_id)
        num_variants = Aligner.get_num_variants(template)

        if previous_num_variants != -1 and \
            num_variants != previous_num_variants:

            raise Exception('Num variants does not match between all \
                templates!')

        previous_num_variants = num_variants


def update_library_alignment_progress(library_progress_string):

    progress_string = alignment_progress_string + '\n' + \
        library_progress_string

    alignment_progress_callback(progress_string)


def set_active_alignment(alignment):
    global active_alignment

    active_alignment = alignment
    db.set_active_alignment(alignment)


def get_active_alignment():

    if active_alignment is None:
        raise Exception('No active alignment specified!')

    return active_alignment


def get_export_path(filename=None):
    if filename is None:
        return get_full_path(export_subdirectory)
    else:
        return os.path.join(get_full_path(export_subdirectory), filename)


def export_csv(filename, header, data):

    export_directory = get_full_path(export_subdirectory)
    mkdir_if_not_exists(export_directory)

    if len(filename) < 4 or filename[-4:] != '.csv':
        filename = filename + '.csv'

    filename = export_directory + "/" + filename

    csv_wrapper.write_csv_file(filename, header, data)


def export_alignment_statistics():

    unique_fields = set()

    # First, get all the headers we'll need
    for alignment in db.get_alignments():

        for library_id, statistics in alignment.statistics.items():

            for statistic_label, statistic in statistics.items():
                unique_fields.add(statistic_label)

    statistic_labels = []

    for statistic_label in unique_fields:
        statistic_labels.append(statistic_label)

    header = ['Alignment', 'Library']
    header.extend(statistic_labels)

    data = []

    for alignment in db.get_alignments():

        for library_id, statistics in alignment.statistics.items():

            row = [alignment.name, db.get_library_by_id(library_id).name]

            for statistic_label in statistic_labels:

                if statistic_label in statistics:
                    row.append(statistics[statistic_label])
                else:
                    row.append('')

            data.append(row)

    export_csv('alignment_statistics', header, data)

# Initialize globals
raw_data_subdirectory = "raw_data"
aligned_subdirectory = ".aligned"
export_subdirectory = "export"
active_alignment = None
alignment_progress_string = ""
fastq_files = {}
experiment_name = None
data_path = None
