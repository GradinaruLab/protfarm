import glob
import os
import json

import Database as db

def get_fastq_files():

    raw_data_directory = get_full_path(raw_data_subdirectory)
    mkdir_if_not_exists(raw_data_directory)

    fastq_files = glob.glob(raw_data_directory + "/*.fastq")

    return fastq_files

def get_alignment_instances():
    pass

def set_workspace_path(new_workspace_path):
    """Set the current workspace path. This is where the db and all other
    files are located"""
    
    global workspace_path

    workspace_path = new_workspace_path

    db.load_database(workspace_path)

def get_full_path(child_path):
    return workspace_path + "/" + child_path

def mkdir_if_not_exists(dir):

    if not os.path.isdir(dir):
        os.mkdir(dir)

# Initialize globals
raw_data_subdirectory = "raw_data"