import Library
from Alignment import Alignment
from Template import Template
import json
import os

def load_database(new_path):

    global path
    global library_db
    global template_db
    global alignment_db
    global library_db_path
    global template_db_path
    global alignment_db_path

    path = new_path

    try:
        library_db_path = path + '/' + library_db_file_name
        template_db_path = path + '/' + template_db_file_name
        alignment_db_path = path + '/' + alignment_db_file_name

        libraries_file = open(library_db_path)
        templates_file = open(template_db_path)
        alignments_file = open(alignment_db_path)

        library_db = json.load(libraries_file)
        template_db = json.load(templates_file)
        alignment_db = json.load(alignments_file)

        print ('Loaded database')
    except IOError:
        initialize_empty_database()
        return

    libraries_file.close()
    templates_file.close()
    alignments_file.close()

def initialize_empty_database():

    global library_db
    global template_db
    global alignment_db

    library_db = {}
    template_db = {}
    alignment_db = {}

    library_db["next_library_id"] = 1
    library_db["libraries"] = {}
    template_db["next_template_id"] = 1
    template_db["templates"] = {}
    alignment_db["next_alignment_id"] = 1
    alignment_db["alignments"] = {}

    dump_database()

def get_libraries():

    library_objects = []

    for library_id, library in library_db["libraries"].items():
        library_object = get_library_object(library_id, library)
        library_objects.append(library_object)

    return library_objects

def get_library(fastq_file):

    for library_id, library in library_db["libraries"].items():
        if fastq_file in library["fastq_files"]:
            return get_library_object(library_id, library)

    raise Exception('FASTQ file doesn\'t exist in any library!')

def add_library(new_library):

    global library_db

    for library_id, library in library_db["libraries"].items():
        if new_library.name == library["name"]:
            raise Exception('Name already exists')

    next_library_id = library_db["next_library_id"]
    library_db["next_library_id"] = next_library_id + 1

    library_db["libraries"][next_library_id] = {}
    library_db["libraries"][next_library_id]["name"] = new_library.name
    library_db["libraries"][next_library_id]["fastq_files"] = []

    update_libraries()

    new_library._id = next_library_id

def remove_library(library):

    if library.id not in library_db["libraries"].keys():
        raise Exception('Library doesn\'t exist in database!')

    del library_db["libraries"][library.id]

    update_libraries()

def update_library(library):

    if library.id not in library_db["libraries"].keys():
        add_library(library)
    else:
        library_db["libraries"][library.id]["name"] = library.name
        library_db["libraries"][library.id]["fastq_files"] = \
            library.fastq_files
        update_libraries()

def get_library_object(library_id, library):
    library_object = Library.Library(library["name"], library_id)
    library_object._fastq_files = library["fastq_files"]
    return library_object

def get_templates():

    template_objects = []

    for template_id, template in template_db["templates"].items():
        template_object = get_template_object(template_id, template)
        template_objects.append(template_object)

    return template_objects

def add_template(new_template):

    global template_db

    for template_id, template in template_db["templates"].items():
        if new_template.sequence == template["sequence"]:
            raise Exception('Template already exists!')

    next_template_id = template_db["next_template_id"]
    template_db["next_template_id"] = next_template_id + 1

    template_db["templates"][next_template_id] = {}
    template_db["templates"][next_template_id]["sequence"] = new_template.sequence

    update_templates()

    new_template._id = next_template_id

def get_template_object(template_id, template):
    template_object = Template.Template(template["sequence"], template_id)
    return template_object

def get_alignments():
    pass

def dump_database():

    update_libraries()
    update_templates()
    update_alignments()

def update_libraries():

    libraries_temp_file = open(library_db_path + '.tmp', 'w')

    json.dump(library_db, libraries_temp_file, indent=4)

    try:
        os.remove(library_db_path)
    except OSError:
        pass

    os.rename(library_db_path + '.tmp', library_db_path)

def update_templates():

    templates_temp_file = open(template_db_path + '.tmp', 'w')

    json.dump(template_db, templates_temp_file, indent=4)

    try:
        os.remove(template_db_path)
    except OSError:
        pass

    os.rename(template_db_path + '.tmp', template_db_path)

def update_alignments():

    alignments_temp_file = open(alignment_db_path + '.tmp', 'w')

    json.dump(alignment_db, alignments_temp_file, indent=4)

    try:
        os.remove(alignment_db_path)
    except OSError:
        pass

    os.rename(alignment_db_path + '.tmp', alignment_db_path)

library_db_file_name = ".libraries.json"
template_db_file_name = ".templates.json"
alignment_db_file_name = ".alignments.json"