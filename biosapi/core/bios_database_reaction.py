import math

def is_number(n):
    try:
        float(n)
    except ValueError:
        return False
    return True

def get_coefficient(i):
    coefficient = None
    if 'coefficient' in i:
        coefficient = i['coefficient']
    elif 'stoichiometry' in i:
        coefficient = str(i['stoichiometry'])
    return coefficient

def dict_to_eq_block(d):
    l = []
    lhs = d
    for id in lhs:
        v = math.fabs(lhs[id])
        if not v == 1:
            l.append("{} {}".format(v, id))
        else:
            l.append(id)
    return ' + '.join(l)

class BiosDatabaseReaction():
    
    def __init__(self, json, api=None, uid=None, id=None, database=None):
        self.json_data=json
        self.api = api
        #self.uid = uid
        #self.id = id
        #self.database = database
    
    @property
    def id(self):
        return self.json_data['entry']
    
    @property
    def uid(self):
        return self.json_data['id']
    
    @property
    def database(self):
        return self.json_data['majorLabel']
    
    @property
    def translocation(self):
        return not len(self.lhs.keys() & self.rhs.keys()) == 0
    
    @property
    def basic(self):
        return len(self.lhs.keys() & self.rhs.keys()) == 0
    
    @property
    def name(self):
        return self.json_data['name']
    
    @property
    def stoichiometry(self):
        s = {}
        for o in self.json_data['leftStoichiometry']:
            v = self.json_data['leftStoichiometry'][o]
            s[o] = -1 * math.fabs(v)
        for o in self.json_data['rightStoichiometry']:
            v = self.json_data['rightStoichiometry'][o]
            s[o] =  1 * math.fabs(v)
        return s
    
    @property
    def cstoichiometry(self):
        stoich = {}
        s = self.stoichiometry
        for id in s:
            p = (int(id), 0)
            if not p in stoich:
                stoich[p] = 0
            stoich[p] += s[id]

        return stoich
    
    def get_raw_stoich(self):
        stoich = {}
        basic = True
        for k in self.json_data['left']:
            coefficient = get_coefficient(self.json_data['left'][k])
            cpd_id = self.json_data['left'][k]['cpdEntry']
            if not cpd_id in stoich:
                stoich[cpd_id] = coefficient
            else:
                basic = False
            basic &= is_number(coefficient)
        for k in self.json_data['right']:
            coefficient = get_coefficient(self.json_data['right'][k])
            cpd_id = self.json_data['right'][k]['cpdEntry']
            if not cpd_id in stoich:
                stoich[cpd_id] = coefficient
            else:
                basic = False
            basic &= is_number(coefficient)
        return stoich, basic
    
    @property
    def lhs(self):
        lhs = {}
        for o in self.json_data['leftStoichiometry']:
            v = self.json_data['leftStoichiometry'][o]
            cpd_data = self.get_internal_compound_data(int(o))
            if cpd_data == None:
                lhs[o] = math.fabs(v)
            else:
                lhs[cpd_data['entry']] = math.fabs(v)
            #print(o, v, self.get_internal_compound_data(int(o)))
        return lhs
    
    @property
    def rhs(self):
        rhs = {}
        for o in self.json_data['rightStoichiometry']:
            v = self.json_data['rightStoichiometry'][o]
            cpd_data = self.get_internal_compound_data(int(o))
            if cpd_data == None:
                rhs[o] = math.fabs(v)
            else:
                rhs[cpd_data['entry']] = math.fabs(v)
            #print(o, v, self.get_internal_compound_data(int(o)))
        return rhs
    
    def get_internal_compound_data(self, id):
        if 'left_component' in self.json_data['connectedEntities']:
            for o in self.json_data['connectedEntities']['left_component']:
                for v in o.values():
                    #print(v['id'], id, v['id'] == id, type(v['id']), type(id))
                    if v['id'] == id:
                        return v
        if 'right_component' in self.json_data['connectedEntities']:
            for o in self.json_data['connectedEntities']['right_component']:
                for v in o.values():
                    if v['id'] == id:
                        return v
        return None
    
    #def __repr__(self):
    #    return self.__str__()
    
    def __str__(self):
        l = dict_to_eq_block(self.lhs)
        r = dict_to_eq_block(self.rhs)
        eq = "<=>"
        return '{}: {} {} {}'.format(self.id, l, eq, r)