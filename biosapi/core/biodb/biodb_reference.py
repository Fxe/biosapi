
class BiodbReference:
    
    def __init__(self, biodb_id, database = ''):
        self._id = biodb_id
        
    def id(self):
        return self._id