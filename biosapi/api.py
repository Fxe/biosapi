from enum import Enum
import requests
import math

	#LigandReaction, 
	#BRENDA,Rhea, BiGG, 
	#@Deprecated BiGG2Reaction,
	#BiGGReaction,
	#Seed, Reactome, 
	#MetaCyc, HumanCyc, EcoCyc, YeastCycReaction, AraCycReaction,
	#PlantCycReaction,
	#BioPath, UniPathway, 
	#ModelSeedReaction,
	#MetaNetXReaction,
	#PIR,
	#NOTFOUND,

class Databases(Enum):
    CHEBI = 'chebi'
    LIPIDMAPS = 'lipidmaps'
    
    BIGG_METABOLITE = 'bigg.metabolite'
    BIGG_REACTION = 'bigg.reaction'
    BIGG1_METABOLITE = 'bigg1.metabolite'
    BIGG1_REACTION = 'bigg1.reaction'
    KEGG_COMPOUND = 'kegg.compound'
    KEGG_GLYCAN = 'kegg.glycan'
    KEGG_DRUG = 'kegg.drug'
    KEGG_REACTION = 'kegg.reaction'
    METACYC_COMPOUND = 'metacyc.compound'
    METACYC_REACTION = 'metacyc.reaction'
    MODELSEED_COMPOUND = 'seed.compound'
    MODELSEED_REACTION = 'seed.reaction'
    
DEFAULT_MAPPER = {
    'BiGG' : Databases.BIGG1_METABOLITE,
    'BiGGMetabolite' : Databases.BIGG_METABOLITE,
    'ChEBI' : Databases.CHEBI,
    'LigandCompound' : Databases.KEGG_COMPOUND,
    'LigandDrug' : Databases.KEGG_DRUG,
    'LigandGlycan' : Databases.KEGG_GLYCAN,
    'LipidMAPS' : Databases.LIPIDMAPS,
    'MetaCyc' : Databases.METACYC_COMPOUND,
    'ModelSeed' : Databases.MODELSEED_COMPOUND,
}

DEFAULT_REACTION_MAPPER = {
    'BiGG' : Databases.BIGG1_REACTION,
    'BiGGReaction' : Databases.BIGG_REACTION,
    'LigandReaction' : Databases.KEGG_REACTION,
    'MetaCyc' : Databases.METACYC_REACTION,
    'ModelSeedReaction' : Databases.MODELSEED_REACTION,
}

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

def get_cprop(cpd, cprop):
    if not cprop in cpd['connectedEntities']:
        return None
    k = set(cpd['connectedEntities'][cprop][0].keys()).pop()
    return cpd['connectedEntities'][cprop][0][k]['properties']['key']

def get_crossreferences(cpd):
    res = set()
    if 'connectedEntities' in cpd and 'has_crossreference_to' in cpd['connectedEntities']:
        for o in cpd['connectedEntities']['has_crossreference_to']:
            v = list(o.values())[0]
            res.add("{}@{}".format(v['entry'], v['majorLabel']))
    return res

def get_cprops(cpd, cprop):
    keys = set()
    if cprop in cpd['connectedEntities']:
        for o in cpd['connectedEntities'][cprop]:
            k = set(o.keys()).pop()
            keys.add(o[k]['properties']['key'])
    return keys

def get_mdl_mol(cpd):
    k = set(cpd['connectedEntities']['has_mdl_mol_file'][0].keys()).pop()
    return cpd['connectedEntities']['has_mdl_mol_file'][0][k]['properties']['key']

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

def get_inchi(cpd):
    inchi = None
    inchikey = None
    if cpd['majorLabel'] == 'LigandCompound':
        mol = get_cprop(cpd, 'has_mdl_mol_file')
        if not mol == None:
            inchi = obb_convert(mol, 'mol', 'inchi').strip()
            inchikey = obb_convert(inchi, 'inchi', 'inchikey').strip()
            return inchi, inchikey
    inchi = get_cprop(cpd, 'has_inchi')
    inchikey = get_cprop(cpd, 'has_inchikey')
    if inchi == None and inchikey == None:
        return None, None
    elif inchi == None or inchikey == None:
        raise Exception('inconsistent')
    inchikey_comp = obb_convert(inchi, 'inchi', 'inchikey').strip()
    if not inchikey == inchikey_comp:
        raise Exception('mismatch keys')
    return inchi, inchikey

        
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
    
class BiosUniversalReaction():
    
    def __init__(self, json, api=None):
        self.json_data=json
        self.api = api
        self.database_reactions = None
        
    @property
    def reactions(self):
        if not self.database_reactions == None:
            return self.database_reactions
        
        ids = set()
        for m in self.json_data['members']:
            ids.add(m['id'])
            
        if not self.api == None:
            self.database_reactions = {}
            for r in self.api.get_reactions(ids):
                self.database_reactions[r.uid] = r
        return self.database_reactions
        
class BIOS():
    
    def __init__(self, base_url = None):
        if base_url:
            self.base_url = base_url
        else:
            self.base_url = 'http://192.168.1.21:8080/biosynth-web-biobase'
        self.headers = {'Content-type': 'application/json', 'Accept': 'application/json'}
    
    def get_rxn_database_ids(self, database):
        resp = requests.get(self.base_url + '/api/dsa/rxn/' + database, headers=self.headers)
        if resp.status_code != 200:
            raise ApiError('GET /api/dsa/rxn/{} {}'.format(database, resp.status_code))
        rxn_ids = set()
        for i in resp.json():
            rxn_ids.add(i['id'])
        return rxn_ids

    def get_cpd_database_ids(self, database):
        resp = requests.get(self.base_url + '/api/dsa/cpd/' + database, headers=self.headers)
        if resp.status_code != 200:
            raise ApiError('GET /api/dsa/cpd/{} {}'.format(database, resp.status_code))
        rxn_ids = set()
        for i in resp.json():
            rxn_ids.add(i['id'])
        return rxn_ids
    
    def get_metabolites(self, ids, database='*'):
        items = []
        for chunk in chunks(list(ids), 200):
        #http://localhost:8080/biosynth-web-biobase//api/dsa/*/cpds
            resp = requests.post(self.base_url + '/api/dsa/' + database + '/cpds', headers=self.headers, json=list(chunk))
            if resp.status_code != 200:
                raise ApiError('GET /api/cura/universal/cpd {}'.format(resp.status_code))
            for item in resp.json():
                items.append(item)
        return items
    
    def get_reactions(self, ids, database='*'):
        items = []
        for chunk in chunks(list(ids), 200):
            #http://localhost:8080/biosynth-web-biobase//api/dsa/*/cpds
            request_url = '{}/api/dsa/{}/rxns'.format(self.base_url, database)
            resp = requests.post(request_url, headers=self.headers, json=chunk)
            
            if resp.status_code != 200:
                raise ApiError('POST {} {}'.format(request_url, resp.status_code))
            for item in resp.json():
                items.append(BiosDatabaseReaction(item, self))
        return items
    
    def get_universal_metabolties(self):
        resp = requests.get(self.base_url + '/api/cura/universal/cpd', headers=self.headers)
        
        if resp.status_code != 200:
            raise ApiError('GET /api/cura/universal/cpd {}'.format(resp.status_code))
            
        return resp.json()
    
    def get_universal_reactions(self):
        resp = requests.get(self.base_url + '/api/cura/universal/rxn', headers=self.headers)
        
        if resp.status_code != 200:
            raise ApiError('GET /api/cura/universal/rxn {}'.format(resp.status_code))
            
        return resp.json()
    
    def get_metabolite(self, id, database='*'):
        resp = requests.get(self.base_url + '/api/dsa/' + database + '/cpd/' + id, headers=self.headers)
        if resp.status_code != 200:
            # This means something went wrong.
            raise ApiError('GET /api/dsa/{}/cpd/ {}'.format(database, resp.status_code))
            
        return resp.json()
    
    def get_model_compartments(self, model_id):
        request_url = '{}/api/model/{}/cmp'.format(self.base_url, model_id)
        resp = requests.get(request_url, headers=self.headers)
        data = resp.json()
        return data
    
    def get_model_species(self, model_id):
        request_url = '{}/api/model/{}/spi'.format(self.base_url, model_id)
        resp = requests.get(request_url, headers=self.headers)
        data = resp.json()
        return data
    
    def get_model_reactions(self, model_id):
        request_url = '{}/api/model/{}/rxn'.format(self.base_url, model_id)
        resp = requests.get(request_url, headers=self.headers)
        data = resp.json()
        return data
    
    def get_model_genes(self, model_id):
        request_url = '{}/api/model/{}/gene'.format(self.base_url, model_id)
        resp = requests.get(request_url, headers=self.headers)
        data = resp.json()
        return data
    
    def get_bios_node_edges(self, node_id):
        request_url = '{}/api/neo4j/data/{}/edges'.format(self.base_url, node_id)
        resp = requests.get(request_url, headers=self.headers)
        data = resp.json()
        return data
    
    def get_bios_node(self, node_id):
        request_url = '{}/api/neo4j/data/{}'.format(self.base_url, node_id)
        resp = requests.get(request_url, headers=self.headers)
        data = resp.json()
        return data