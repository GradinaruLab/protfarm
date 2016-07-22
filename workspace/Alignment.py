
class Alignment(object):

    def __init__(self, method, parameters, library_templates, id = 0):
        self_.method = method

        if id == 0:
            db.add_alignment(self)
        else:
            self._id = id

    @property
    def name(self):
        return self._name

    @property
    def method(self):
        return self._method
    
    @property
    def parameters(self):
        return self._parameters
    
    @property
    def library_templates(self):
        return self._library_templates
