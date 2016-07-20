import os

class Library:

    def __init__(self, name):
        self.fastq_files = []
        self.name = name

    def add_fastq_file(self, file_name):
        self.fastq_files.append(file_name)

    def get_fastq_files(self):
        return self.fastq_files
        
    @staticmethod
    def from_file(file_name):

        path, file_tail = os.path.split(file_name)

        library_name = file_tail[:-4]

        library = Library(library_name)

        return library