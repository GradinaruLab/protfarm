import Database as db

class Template(object):

    def __init__(self, sequence, id = 0):
        self._sequence = sequence

        if id == 0:
            db.add_template(self)
        else:
            self._id = id

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