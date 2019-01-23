
class FASTQ_File(object):

    def __init__(self, name, id = None):

        from . import Database as db

        self._name = name
        self._reverse_complement = False

        if id == None:
            db.add_FASTQ_file(self)
        else:
            self._id = id

    @property
    def is_reverse_complement(self):
        is_reverse_complement = self._reverse_complement

        return is_reverse_complement

    @property
    def name(self):
        return self._name

    @property
    def id(self):
        return self._id

    @is_reverse_complement.setter
    def is_reverse_complement(self, reverse_complement):

        from . import Database as db
        self._reverse_complement = reverse_complement
        db.update_FASTQ_file(self)

    @id.setter
    def id(self, new_id):
        raise Exception('Can\'t change id of a FASTQ File!')

    @name.setter
    def name(self, new_id):
        raise Exception('Can\'t change name of a FASTQ File!')