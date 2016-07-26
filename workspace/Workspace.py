import glob
import os
import json

import Database as db

from fileio import csv_wrapper

def get_fastq_files():

    raw_data_directory = get_full_path(raw_data_subdirectory)
    mkdir_if_not_exists(raw_data_directory)

    fastq_files = glob.glob(raw_data_directory + "/*.fastq")

    return fastq_files

def write_sequence_file(library, alignment, sequence_uuids):

    file_name = workspace_path + "/" + aligned_subdirectory + "/" + \
        str(library.id) + "_" + str(alignment.id) + ".csv"

    csv_wrapper.write_csv_file(file_name, ['Sequence', 'UUID'], sequence_uuids)

def get_alignment_instances():
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

def get_raw_data_path(child_path):
    return workspace_path + "/" + raw_data_subdirectory + "/" + child_path

def get_full_path(child_path):
    return workspace_path + "/" + child_path

def mkdir_if_not_exists(dir):

    if not os.path.isdir(dir):
        os.mkdir(dir)

# Initialize globals
raw_data_subdirectory = "raw_data"
aligned_subdirectory = ".aligned"