

class BiosModelReaction:
    
    def __init__(self, data):
        self.data = data
    
    @property
    def id(self):
        return self.data['id']
    
    @property
    def uid(self):
        return self.data['bios_id']
    
    @property
    def name(self):
        return self.data['name']
    
    def decode_stoich_value(self, v):
        if v == None or len(v.strip()) == 0:
            return 1
        else:
            return float(v)
    
    def get_stoichiometry(self):
        s = {}
        for o in self.data['bios_stoichiometry']['l']:
            s[o[0]] = -1 * self.decode_stoich_value(o[2])
        for o in self.data['bios_stoichiometry']['r']:
            s[o[0]] = self.decode_stoich_value(o[2])
        return s

    def get_cstoichiometry(self, metabolites):
        s = {}
        for o in self.data['bios_stoichiometry']['l']:
            if o[0] in metabolites:
                s[(o[0], metabolites[o[0]]['compartment'])] = -1 * self.decode_stoich_value(o[2])
        for o in self.data['bios_stoichiometry']['r']:
            if o[0] in metabolites:
                s[(o[0], metabolites[o[0]]['compartment'])] = self.decode_stoich_value(o[2])
        return s
