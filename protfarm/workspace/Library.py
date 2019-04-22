import os
from . import FASTQ_File

class Library(object):

    def __init__(self, name, id=0):

        from . import Database as db

        self._fastq_files = []
        self._name = name
        self._metadata = {}

        if id == 0:
            db.add_library(self)
        else:
            self._id = id

    def add_file(self, file_name):

        file_name = FASTQ_File.clean_FASTQ_file_name(file_name)

        from . import Database as db

        try:
            existing_file_index = self._fastq_files.index(file_name)
            return
        except:
            pass

        self._fastq_files.append(file_name)
        db.update_library(self)

    def remove_file(self, file_name):

        from . import Database as db
        from . import Workspace as ws

        file_to_remove_index = self._fastq_files.index(file_name)
        del self._fastq_files[file_to_remove_index]
        db.update_library(self)
        ws.remove_library_alignments(self)

    @property
    def fastq_files(self):
        copy = self._fastq_files[:]
        return copy

    @fastq_files.setter
    def fastq_files(self, new_fastq_files):

        for new_fastq_file in new_fastq_files:
            self.add_file(new_fastq_file)

        files_to_remove = []

        for existing_fastq_file in self._fastq_files:
            if existing_fastq_file not in new_fastq_files:
                files_to_remove.append(existing_fastq_file)

        for existing_fastq_file in files_to_remove:
            self.remove_file(existing_fastq_file)

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, new_name):

        from . import Database as db
        self._name = new_name
        db.update_library(self)

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self, new_id):
        raise Exception('Can\'t change id of a library!')

    @property
    def metadata(self):
        return self._metadata

    @metadata.setter
    def metadata(self, new_metadata):

        from . import Database as db
        self._metadata = new_metadata
        db.update_library(self)
