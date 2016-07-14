import glob
import os

from Library import Library

raw_data_subdirectory = "raw_data"
library_subdirectory = ".libraries"
workspace_path = "."

def get_fastq_files():

    raw_data_directory = get_full_path(raw_data_subdirectory)
    mkdir_if_not_exists(raw_data_directory)

    fastq_files = glob.glob(raw_data_directory + "/*.fastq")

    return fastq_files

def get_libraries():

    library_directory = get_full_path(library_subdirectory)
    mkdir_if_not_exists(library_directory)

    library_files = glob.glob(library_directory + "/*.lib")

    libraries = []

    for library_file in library_files:
        libraries.append(Library.from_file(library_file))

    return libraries

def set_workspace_path(path):
    global workspace_path
    workspace_path = path

def get_full_path(subdirectory):
    return workspace_path + "/" + subdirectory

def mkdir_if_not_exists(dir):

    if not os.path.isdir(dir):
        os.mkdir(dir)