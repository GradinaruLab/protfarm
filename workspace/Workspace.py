import glob
import os
import json
import importlib
import subprocess
from . import Database as db
from .FASTQ_File import FASTQ_File

from fileio import csv_wrapper

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

    if os.path.isfile(fastq_file_path):
        fastq_file = open(fastq_file_path, 'r')
    elif os.path.isfile(fastq_file_path + '.gz'):
        gunzip_command = 'gzip -dk ' + fastq_file_path + '.gz'
        subprocess.call(gunzip_command, shell=True)
        fastq_file = open(fastq_file_path, 'r')
    else:
        raise Exception("File does not exist!")

    fastq_files[fastq_file_name] = fastq_file

    return fastq_file

def close_fastq_file(fastq_file_name):

    fastq_files[fastq_file_name].close()

    fastq_file_path = get_raw_data_path(fastq_file_name)

    if os.path.isfile(fastq_file_path + '.gz'):
        os.remove(fastq_file_path)

def get_full_path(child_path):
    return os.path.join(workspace_path, child_path)

def mkdir_if_not_exists(dir):

    if not os.path.isdir(dir):
        os.makedirs(dir)

def align_all(callback):

    from sequencing.Perfect_Match_Aligner import Perfect_Match_Aligner
    from sequencing.Bowtie_Aligner import Bowtie_Aligner
    from sequencing.Aligner import Aligner

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

        Aligner.validate_alignment(alignment)

        for library_id, template_id in alignment.library_templates.items():

            library = db.get_library_by_id(library_id)

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

            aligner.align(alignment, library, \
                update_library_alignment_progress)

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

def get_export_path():
    return get_full_path(export_subdirectory)

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
