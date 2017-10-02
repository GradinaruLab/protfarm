from utils import DNA

class Template(object):

    def __init__(self, sequence, name = None, id = 0, \
        reverse_complement_template_id = None):

        from . import Database as db
        self._sequence = sequence
        self._reverse_complement_template_id = None

        if name == None:
            self._name = sequence
        else:
            self._name = name

        if id == 0:
            db.add_template(self)
            self.reverse_complement_template_id = reverse_complement_template_id
        else:
            self._id = id
            self._reverse_complement_template_id = reverse_complement_template_id

    def get_variant_template(self):

        variant_template = ''

        for sequence_element in self._sequence:
            if sequence_element in DNA.IUPAC and \
                sequence_element not in DNA.get_nucleotides():
                variant_template += sequence_element

        return variant_template

    @property
    def name(self):
        return self._name

    @property
    def reverse_complement_template_id(self):
        return self._reverse_complement_template_id

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, new_sequence):
        if new_sequence == self._sequence:
            return

        raise Exception('Can\'t change sequence of a template!')

    @property
    def id(self):
        return self._id

    @name.setter
    def name(self, new_name):
        if new_name == self._name:
            return

        self._name = new_name

        from . import Database as db
        db.update_template(self)

    @id.setter
    def id(self, new_id):
        if new_id == self._id:
            return

        raise Exception('Can\'t change id of a template!')

    @reverse_complement_template_id.setter
    def reverse_complement_template_id(self, reverse_complement_template_id):

        if self._reverse_complement_template_id == reverse_complement_template_id:
            return
        
        from . import Database as db

        # If the previous reverse complement template wasn't empty, let's update it
        if self._reverse_complement_template_id != None:
            other_template = db.get_template_by_id(int(self._reverse_complement_template_id))
            other_template._reverse_complement_template_id = None
            db.update_template(other_template)

        self._reverse_complement_template_id = reverse_complement_template_id
        db.update_template(self)

        # If the new template isn't empty, let's update it
        if self._reverse_complement_template_id != None:
            other_template = db.get_template_by_id(int(reverse_complement_template_id))
            other_template._reverse_complement_template_id = self._id
            db.update_template(other_template)