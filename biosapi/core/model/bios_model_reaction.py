class BiosModelReaction():
    
    def __init__(self, json, api=None, uid=None, id=None, model=None):
        self.json_data = json
        self.api = api
        self.uid = self.json_data['bios_id']
        self.id = self.json_data['id']
        #self.database = database
        
    @property
    def cstoichiometry(self):
        cstoichiometry = {}
        if 'bios_stoichiometry' in self.json_data:
            for h in self.json_data['bios_stoichiometry']['l']:
                cstoichiometry[(h[0], '?')] = -1 * float(h[2])
            for h in self.json_data['bios_stoichiometry']['r']:
                cstoichiometry[(h[0], '?')] = float(h[2])
        return cstoichiometry