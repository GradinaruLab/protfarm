import glob
import os
import json
import importlib

import Database as db
from fileio import csv_wrapper

def get_fastq_files():

    raw_data_directory = get_full_path(raw_data_subdirectory)
    mkdir_if_not_exists(raw_data_directory)

    # fastq_files = glob.glob(raw_data_directory + "/*.fastq")
    fastq_files = [file for file in os.listdir(raw_data_directory) \
        if file.endswith('.fastq')]

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

def set_workspace_path(new_workspace_path):
    """Set the current workspace path. This is where the db and all other
    files are located"""
    
    global workspace_path

    workspace_path = new_workspace_path

    try:
        os.mkdir(raw_data_subdirectory)
    except:
        pass

    try:
        os.mkdir(aligned_subdirectory)
    except:
        pass

    db.load_database(workspace_path)

    global active_alignment

    try:
        active_alignment = db.get_active_alignment()
    except:
        active_alignment = None

def get_raw_data_path(child_path):
    return workspace_path + "/" + raw_data_subdirectory + "/" + child_path

def get_full_path(child_path):
    return workspace_path + "/" + child_path

def mkdir_if_not_exists(dir):

    if not os.path.isdir(dir):
        os.mkdir(dir)

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
        if alignment.method not in [Perfect_Match_Aligner.__name__]:
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

# Initialize globals
raw_data_subdirectory = "raw_data"
aligned_subdirectory = ".aligned"
export_subdirectory = "export"
active_alignment = None
alignment_progress_string = ""