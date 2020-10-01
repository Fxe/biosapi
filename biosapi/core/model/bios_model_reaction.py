

class BiosModelReaction:
    
    def __init__(self, json, api=None, uid=None, id=None, model=None):
        self.json_data = json
        self.api = api
        #self.database = database

    @property
    def id(self):
        return self.json_data['id']

    @property
    def uid(self):
        return self.json_data['bios_id']

    @property
    def name(self):
        return self.json_data['name']

    def decode_stoich_value(self, v):
        if v == None or len(v.strip()) == 0:
            return 1
        else:
            return float(v)

    def get_stoichiometry(self):
        s = {}
        for o in self.json_data['bios_stoichiometry']['l']:
            s[o[0]] = -1 * self.decode_stoich_value(o[2])
        for o in self.json_data['bios_stoichiometry']['r']:
            s[o[0]] = self.decode_stoich_value(o[2])
        return s

    def get_cstoichiometry(self, metabolites):
        s = {}
        for o in self.json_data['bios_stoichiometry']['l']:
            if o[0] in metabolites:
                s[(o[0], metabolites[o[0]]['compartment'])] = -1 * self.decode_stoich_value(o[2])
        for o in self.json_data['bios_stoichiometry']['r']:
            if o[0] in metabolites:
                s[(o[0], metabolites[o[0]]['compartment'])] = self.decode_stoich_value(o[2])
        return s
        
    @property
    def cstoichiometry(self):
        cstoichiometry = {}
        if 'bios_stoichiometry' in self.json_data:
            for h in self.json_data['bios_stoichiometry']['l']:
                if type(h[2]) == str and len(h[2]) == 0:
                    h[2] = '1'
                cstoichiometry[(h[0], '?')] = -1 * float(h[2])
            for h in self.json_data['bios_stoichiometry']['r']:
                if type(h[2]) == str and len(h[2]) == 0:
                    h[2] = '1'
                cstoichiometry[(h[0], '?')] = float(h[2])
        return cstoichiometry
