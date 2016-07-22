import os
import Database as db

class Library(object):

    def __init__(self, name, id = 0):
        self._fastq_files = []
        self._name = name

        if id == 0:
            db.add_library(self)
        else:
            self._id = id

    def append_file(self, file_name):
        self._fastq_files.append(file_name)
        db.update_library(self)

    def extend_files(self, file_names):
        self._fastq_files.extend(file_names)
        db.update_library(self)

    @property
    def fastq_files(self):
        copy = self._fastq_files[:]
        return copy

    @fastq_files.setter
    def fastq_files(self, new_fastq_files):
        self._fastq_files = new_fastq_files
        db.update_library(self)

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, new_name):
        self._name = new_name
        db.update_library(self)

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self, new_id):
        raise Exception('Can\'t change id of a library!')
    