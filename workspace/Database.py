from . import Library
from . import Alignment
from . import Template
import json
import os

from utils import utils

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
    alignment_db["active_alignment"] = 0

    dump_database()

def get_libraries():

    library_objects = []

    for library_id, library in library_db["libraries"].items():
        library_object = get_library_object(library_id, library)
        library_objects.append(library_object)

    return library_objects

def get_associated_library(fastq_file):

    for library_id, library in library_db["libraries"].items():
        if fastq_file in library["fastq_files"]:
            return get_library_object(library_id, library)

    raise Exception('FASTQ file doesn\'t exist in any library!')

def get_library(name):

    for library_id, library in library_db["libraries"].items():
        if name == library["name"]:
            return get_library_object(library_id, library)

    raise Exception('No library found with name \'' + name + '\'')

def get_library_by_id(id):

    if str(id) not in library_db["libraries"].keys():
        raise Exception('No library with id \'' + str(id) + '\'')

    return get_library_object(id, library_db["libraries"][str(id)])

def add_library(new_library):

    global library_db

    for library_id, library in library_db["libraries"].items():
        if new_library.name == library["name"]:
            raise Exception('Name already exists')

    next_library_id = library_db["next_library_id"]
    library_db["next_library_id"] = next_library_id + 1

    library_db["libraries"][str(next_library_id)] = {}
    library_db["libraries"][str(next_library_id)]["name"] = new_library.name
    library_db["libraries"][str(next_library_id)]["fastq_files"] = []

    update_libraries()

    new_library._id = next_library_id

def remove_library(library):

    if str(library.id) not in library_db["libraries"].keys():
        raise Exception('Library doesn\'t exist in database!')

    del library_db["libraries"][str(library.id)]

    update_libraries()

def update_library(library):

    if str(library.id) not in library_db["libraries"].keys():
        add_library(library)
    else:
        library_db["libraries"][str(library.id)]["name"] = library.name
        library_db["libraries"][str(library.id)]["fastq_files"] = \
            library.fastq_files
        update_libraries()

def get_library_object(library_id, library):
    library_object = Library.Library(library["name"], int(library_id))
    library_object._fastq_files = library["fastq_files"]
    return library_object

def get_template_seed():
    return template_db["next_template_id"]

def get_library_seed():
    return library_db["next_library_id"]

def get_alignment_seed():
    return alignment_db["next_alignment_id"]

def get_templates():

    template_objects = []

    for template_id, template in template_db["templates"].items():
        template_object = get_template_object(template_id, template)
        template_objects.append(template_object)

    return template_objects

def get_template_by_sequence(sequence):
    for template_id, template in template_db['templates'].items():
        if sequence == template["sequence"]:
            return get_template_object(template_id, template)
    return None

def get_template_by_id(id):

    if str(id) not in template_db["templates"].keys():
        raise Exception('No template with id \'' + str(id) + '\'')

    return get_template_object(id, template_db["templates"][str(id)])

def add_template(new_template):

    global template_db

    for template_id, template in template_db["templates"].items():
        if new_template.sequence == template["sequence"]:
            raise Exception('Template already exists!')

    next_template_id = template_db["next_template_id"]
    template_db["next_template_id"] = next_template_id + 1

    template_db["templates"][str(next_template_id)] = {}
    template_db["templates"][str(next_template_id)]["sequence"] = \
        new_template.sequence

    update_templates()

    new_template._id = next_template_id

def get_template_object(template_id, template):

    template_object = Template.Template(template["sequence"], int(template_id))
    return template_object

def get_alignments():

    alignment_objects = []

    for alignment_id, alignment in alignment_db["alignments"].items():
        alignment_object = get_alignment_object(alignment_id, alignment)
        alignment_objects.append(alignment_object)

    return alignment_objects

def get_alignment(name):

    for alignment in get_alignments():
        if alignment.name == name:
            return alignment

    raise Exception('No alignment exists with name \'' + name + '\'')

def get_alignment_by_id(id):

    if str(id) not in alignment_db["alignments"].keys():
        raise Exception('No alignment with id \'' + str(id) + '\'')

    return get_alignment_object(id, alignment_db["alignments"][str(id)])

def add_alignment(new_alignment):

    global alignment_db

    next_alignment_id = alignment_db["next_alignment_id"]

    alignment_db["next_alignment_id"] = next_alignment_id + 1

    alignment_db["alignments"][str(next_alignment_id)] = {}
    alignment_db["alignments"][str(next_alignment_id)]["method"] = \
        new_alignment.method
    alignment_db["alignments"][str(next_alignment_id)]["parameters"] = \
        new_alignment.parameters
    alignment_db["alignments"][str(next_alignment_id)]["library_templates"] = \
        new_alignment.library_templates
    alignment_db["alignments"][str(next_alignment_id)]["statistics"] = \
        new_alignment.statistics

    update_alignments()

    new_alignment._id = next_alignment_id

def get_alignment_object(alignment_id, alignment):

    library_templates = \
        utils.convert_string_keys_to_ints(alignment["library_templates"])

    statistics = \
        utils.convert_string_keys_to_ints(alignment["statistics"])

    alignment_object = Alignment.Alignment(alignment["method"], \
        alignment["parameters"], library_templates, statistics, \
        int(alignment_id))

    return alignment_object

def update_alignment(alignment):

    if str(alignment.id) not in alignment_db["alignments"].keys():
        add_alignment(alignment)
    else:
        alignment_db["alignments"][str(alignment.id)]["method"] = \
            alignment.method
        alignment_db["alignments"][str(alignment.id)]["parameters"] = \
            alignment.parameters
        alignment_db["alignments"][str(alignment.id)]["library_templates"] = \
            alignment.library_templates
        alignment_db["alignments"][str(alignment.id)]["statistics"] = \
            alignment.statistics

        update_alignments()

def set_active_alignment(alignment):

    alignment_db["active_alignment"] = alignment.id

    update_alignments()

def get_active_alignment():

    return get_alignment_by_id(alignment_db["active_alignment"])

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
