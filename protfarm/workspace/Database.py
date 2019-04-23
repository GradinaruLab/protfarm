import json
import os
import shutil

from pepars.utils import utils

from . import Library
from . import Alignment
from . import Template
from . import FASTQ_File


def load_database(new_path):

    global path
    global metadata
    global library_db
    global template_db
    global alignment_db
    global metadata_path
    global library_db_path
    global template_db_path
    global alignment_db_path

    path = new_path

    try:
        metadata_path = os.path.join(path, metadata_file_name)
        library_db_path = os.path.join(path, library_db_file_name)
        template_db_path = os.path.join(path, template_db_file_name)
        alignment_db_path = os.path.join(path, alignment_db_file_name)

        if not os.path.exists(metadata_path):
            initialize_metadata_file()

        metadata_file = open(metadata_path)
        libraries_file = open(library_db_path)
        templates_file = open(template_db_path)
        alignments_file = open(alignment_db_path)

        metadata = json.load(metadata_file)
        library_db = json.load(libraries_file)
        template_db = json.load(templates_file)
        alignment_db = json.load(alignments_file)
		
        metadata_file.close()
        libraries_file.close()
        templates_file.close()
        alignments_file.close()

        update_database()

    except IOError:
        initialize_empty_database()
        return

    try:
        metadata_file.close()
    except Exception:
        pass
    try:
        libraries_file.close()
    except Exception:
        pass

    try:
        templates_file.close()
    except Exception:
        pass
    
    try:
        alignments_file.close()
    except Exception:
        pass

def initialize_metadata_file():

    global metadata

    metadata = {}

    metadata["version"] = 0

    update_metadata()

def initialize_empty_database():

    global metadata
    global library_db
    global template_db
    global alignment_db

    metadata = {}
    library_db = {}
    template_db = {}
    alignment_db = {}

    metadata["version"] = LATEST_VERSION
    library_db["next_library_id"] = 1
    library_db["libraries"] = {}
    library_db["next_FASTQ_file_id"] = 1
    library_db["FASTQ_files"] = {}
    template_db["next_template_id"] = 1
    template_db["templates"] = {}
    alignment_db["next_alignment_id"] = 1
    alignment_db["alignments"] = {}
    alignment_db["active_alignment"] = 0

    dump_database()

def get_samples():
    return get_libraries()


def get_sample_names():
    return [x.name for x in get_libraries()]


def get_libraries(metadata_filter=None):

    if metadata_filter is None:
        metadata_filter = {}

    library_objects = []

    for library_id, library in library_db["libraries"].items():
        library_object = get_library_object(library_id, library)

        passes_filter = True

        for key, value in metadata_filter.items():
            if key not in library_object.metadata:
                passes_filter = False
                break
            if library_object.metadata[key] != value:
                passes_filter = False
                break

        if passes_filter:
            library_objects.append(library_object)

    library_objects = sorted(library_objects, key=lambda x: x.name)

    return library_objects

def get_FASTQ_files():

    FASTQ_file_objects = []

    if "FASTQ_files" not in library_db:
        return []

    for FASTQ_file_id, FASTQ_file in sorted(library_db["FASTQ_files"].items(), key=lambda x: x[1]["name"]):
        FASTQ_file_object = get_FASTQ_file_object(FASTQ_file_id, FASTQ_file)
        FASTQ_file_objects.append(FASTQ_file_object)

    return FASTQ_file_objects

def get_FASTQ_file(name):

    for FASTQ_file_id, FASTQ_file in library_db["FASTQ_files"].items():
        if name == FASTQ_file["name"]:
            return get_FASTQ_file_object(FASTQ_file_id, FASTQ_file)

    raise Exception('No FASTQ file found with name \'' + name + '\'')

def get_associated_library(fastq_file):

    for library_id, library in library_db["libraries"].items():
        if fastq_file in library["fastq_files"]:
            return get_library_object(library_id, library)

    return None

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
    library_db["libraries"][str(next_library_id)]["metadata"] = {}

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
        library_db["libraries"][str(library.id)]["metadata"] = library.metadata
        update_libraries()

def update_FASTQ_file(FASTQ_file):

    if str(FASTQ_file.id) not in library_db["FASTQ_files"].keys():
        add_FASTQ_file(FASTQ_file)
    else:
        library_db["FASTQ_files"][str(FASTQ_file.id)]["name"] = FASTQ_file.name
        library_db["FASTQ_files"][str(FASTQ_file.id)]["reverse_complement"] = \
            FASTQ_file.is_reverse_complement
        update_libraries()

def get_library_object(library_id, library):
    library_object = Library(library["name"], int(library_id))
    library_object._fastq_files = library["fastq_files"]
    if "metadata" in library:
        library_object._metadata = library["metadata"]
    return library_object

def get_FASTQ_file_object(FASTQ_file_id, FASTQ_file):
    FASTQ_file_object = FASTQ_File(FASTQ_file["name"], int(FASTQ_file_id))
    FASTQ_file_object.is_reverse_complement = FASTQ_file["reverse_complement"]
    return FASTQ_file_object

def get_template_seed():
    return template_db["next_template_id"]

def get_library_seed():
    return library_db["next_library_id"]

def get_alignment_seed():
    return alignment_db["next_alignment_id"]

def get_templates():

    template_objects = []

    for template_id, template in sorted(template_db["templates"].items()):
        template_object = get_template_object(template_id, template)
        template_objects.append(template_object)

    return template_objects

def get_template_by_sequence(sequence):
    for template_id, template in template_db['templates'].items():
        if sequence == template["sequence"]:
            return get_template_object(template_id, template)
    return None

def get_template_by_name(name):
    for template_id, template in template_db['templates'].items():
        if name == template["name"]:
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
        if new_template.name == template["name"]:
            raise Exception("Template with name '%s' already exists!" % 
                new_template.name)

    next_template_id = template_db["next_template_id"]
    template_db["next_template_id"] = next_template_id + 1

    template_db["templates"][str(next_template_id)] = {}
    template_db["templates"][str(next_template_id)]["sequence"] = \
        new_template.sequence
    template_db["templates"][str(next_template_id)]["reverse_complement_template_id"] = new_template.reverse_complement_template_id
    template_db["templates"][str(next_template_id)]["name"] = new_template.name

    update_templates()

    new_template._id = next_template_id

def delete_template(template):
    # Get its reverse complement to set it to null also

    if template.reverse_complement_template_id != None:
        reverse_complement_template = get_template_by_id(\
            template.reverse_complement_template_id)
        reverse_complement_template.reverse_complement_template_id = None

    del template_db["templates"][str(template.id)]

    update_templates()

def delete_sample(sample):

    del library_db["libraries"][str(sample.id)]

    update_libraries()

def update_template(template):

    if str(template.id) not in template_db["templates"].keys():
        add_template(template)
    else:
        template_db["templates"][str(template.id)]["sequence"] = template.sequence
        template_db["templates"][str(template.id)]["reverse_complement_template_id"] = \
            template.reverse_complement_template_id
        template_db["templates"][str(template.id)]["name"] = template.name
        update_templates()

def get_template_object(template_id, template):

    if "reverse_complement_template_id" in template:
        reverse_complement_template_id = template["reverse_complement_template_id"]
    else:
        reverse_complement_template_id = None

    template_object = Template(template["sequence"], id=int(template_id),\
        name=template["name"], reverse_complement_template_id = reverse_complement_template_id)

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

def add_FASTQ_file(new_FASTQ_file):

    global library_db

    if "next_FASTQ_file_id" not in library_db:
        next_FASTQ_file_id = 1
    else:
        next_FASTQ_file_id = library_db["next_FASTQ_file_id"]

    if "FASTQ_files" not in library_db:
        library_db["FASTQ_files"] = {}

    library_db["next_FASTQ_file_id"] = next_FASTQ_file_id + 1

    for existing_FASTQ_file in library_db["FASTQ_files"].values():
        if existing_FASTQ_file["name"] == new_FASTQ_file.name:
            raise Exception("Can't add FASTQ file '%s' because existing FASTQ file with that name exists!" % new_FASTQ_file.name)

    library_db["FASTQ_files"][str(next_FASTQ_file_id)] = {}
    library_db["FASTQ_files"][str(next_FASTQ_file_id)]["name"] = new_FASTQ_file.name
    library_db["FASTQ_files"][str(next_FASTQ_file_id)]["reverse_complement"] = new_FASTQ_file.is_reverse_complement

    update_libraries()

    new_FASTQ_file._id = next_FASTQ_file_id


def add_alignment(new_alignment):

    global alignment_db

    alignment = None

    try:
        alignment = get_alignment_by_parameters(
            new_alignment.method,
            new_alignment.parameters,
            new_alignment.library_templates)
    except ValueError:
        pass

    if alignment:
        raise ValueError("Trying to create alignment that already exists!")

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


def get_alignment_by_parameters(method, parameters, library_templates):

    for alignment in get_alignments():

        if alignment.method == method and \
                alignment.parameters == parameters and \
                alignment.library_templates == library_templates:
            return alignment

    raise ValueError('No alignment exists with these parameters')


def get_alignment_object(alignment_id, alignment):

    library_templates = \
        utils.convert_string_keys_to_ints(alignment["library_templates"])

    statistics = \
        utils.convert_string_keys_to_ints(alignment["statistics"])

    alignment_object = Alignment(alignment["method"], \
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

    update_metadata()
    update_libraries()
    update_templates()
    update_alignments()

def update_metadata():

    metadata_temp_file = open(metadata_path + ".tmp", "w")

    json.dump(metadata, metadata_temp_file, indent=4)
	
    metadata_temp_file.close()

    try:
        os.remove(metadata_path)
    except OSError:
        pass

    shutil.move(metadata_path + ".tmp", metadata_path)

def update_libraries():

    libraries_temp_file = open(library_db_path + '.tmp', 'w')

    json.dump(library_db, libraries_temp_file, indent=4)
	
    libraries_temp_file.close()

    try:
        os.remove(library_db_path)
    except OSError:
        pass

    shutil.move(library_db_path + '.tmp', library_db_path)

def update_templates():

    templates_temp_file = open(template_db_path + '.tmp', 'w')

    json.dump(template_db, templates_temp_file, indent=4)
	
    templates_temp_file.close()

    try:
        os.remove(template_db_path)
    except OSError:
        pass

    shutil.move(template_db_path + '.tmp', template_db_path)

def update_alignments():

    alignments_temp_file = open(alignment_db_path + '.tmp', 'w')

    json.dump(alignment_db, alignments_temp_file, indent=4)
	
    alignments_temp_file.close()

    try:
        os.remove(alignment_db_path)
    except OSError:
        pass

    shutil.move(alignment_db_path + '.tmp', alignment_db_path)

def get_database_version():

    return metadata["version"]

def update_database():

    version = get_database_version()

    while version < LATEST_VERSION:

        update_database_from_version(version)
        version = get_database_version()


def update_database_from_version(current_version):

    if current_version == 0:
        add_names_to_templates()

    metadata["version"] = current_version + 1
    update_metadata()

def add_names_to_templates():

    for template_id, template in sorted(template_db["templates"].items()):
        template_db["templates"][template_id]["name"] = str(template_id)

    update_templates()

LATEST_VERSION = 1
metadata_file_name = ".db.json"
library_db_file_name = ".libraries.json"
template_db_file_name = ".templates.json"
alignment_db_file_name = ".alignments.json"
