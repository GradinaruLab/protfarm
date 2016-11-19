from utils import DNA

class Template(object):

    def __init__(self, sequence, id = 0):

        from . import Database as db
        self._sequence = sequence

        if id == 0:
            db.add_template(self)
        else:
            self._id = id

    def get_variant_template(self):

        variant_template = ''

        for sequence_element in self._sequence:
            if sequence_element in DNA.IUPAC and \
                sequence_element not in DNA.get_nucleotides():
                variant_template += sequence_element

        return variant_template

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, new_sequence):
        raise Exception('Can\'t change sequence of a template!')

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self, new_id):
        raise Exception('Can\'t change id of a template!')