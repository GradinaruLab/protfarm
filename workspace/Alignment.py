import Database as db

class Alignment(object):

    def __init__(self, method, parameters, library_templates, \
        statistics = {}, id = 0):

        self._method = method
        self._parameters = parameters
        self._library_templates = library_templates
        self._statistics = statistics

        if id == 0:
            db.add_alignment(self)
        else:
            self._id = id

    def add_statistics(self, library, statistics):
        self._statistics[library.id] = statistics
        db.update_alignment(self)

    def remove_library(self, library):
        if library.id in self._statistics:
            del self._statistics[library.id]
            db.update_alignment(self)

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

    @property
    def statistics(self):
        return self._statistics

    @property
    def id(self):
        return self._id
    
    @method.setter
    def method(self, new_method):
        raise Exception('You can\'t do that!')

    @parameters.setter
    def parameters(self, new_parameters):
        raise Exception('You can\'t do that!')

    @library_templates.setter
    def library_templates(self, new_library_templates):
        raise Exception('You can\'t do that!')