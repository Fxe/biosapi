import logging, json
from cobra import Metabolite, Reaction, Model

logger = logging.getLogger(__name__)


class BiosModelToCobraBuilder:
    
    def __init__(self, model_cmps, model_spis, model_rxns, model_genes):
        self.model_cmps = model_cmps
        self.model_spis = model_spis
        self.model_rxns = model_rxns
        self.model_genes = model_genes

    @staticmethod
    def from_api(model_id, api):
        model_cmps = api.get_model_compartments(model_id)
        model_spis = api.get_model_species(model_id)
        model_rxns = api.get_model_reactions(model_id)
        model_genes = api.get_model_genes(model_id)
        return BiosModelToCobraBuilder(model_cmps, model_spis, model_rxns, model_genes)

    def convert_modelcompound(self, m):
        mc_id = m['id'] if 'id' in m else "bios_{}".format(m['bios_id'])
        name = m['name'] if 'name' in m else ""
        formula = m['chemicalFormula'] if 'chemicalFormula' in m else ''
        #charge = get_int('charge', 0, metabolite.data)
        #mc_id = metabolite.id
        #annotation = {}
        #if 'dblinks' in metabolite.data:
        #    annotation = get_cpd_annotation(metabolite.data['dblinks'])
        compartment = m['compartment']
        id = mc_id

        met = Metabolite(id, 
                         formula=formula, 
                         name=name, 
                         charge=0, 
                         compartment=compartment)
        #met.annotation[SBO_ANNOTATION] = "SBO:0000247" #simple chemical - Simple, non-repetitive chemical entity.
        #if id.startswith('cpd'):
        #    met.annotation["seed.compound"] = id.split("_")[0]
        #met.annotation.update(annotation)
        return met
    
    def convert_modelreaction_stoichiometry(self, reaction):
        object_stoichiometry = {}
        s = reaction['bios_stoichiometry']
        for metabolite_id, bios_id, v in s['l']:
            if not v:
                v = 1
            if bios_id in self.metabolites:
                object_stoichiometry[self.metabolites[bios_id]] = -1 * float(v)
                #print(metabolites[metabolite_id])
            else:
                print('!')
            pass
        for metabolite_id, bios_id, v in s['r']:
            if not v:
                v = 1
            if bios_id in self.metabolites:
                object_stoichiometry[self.metabolites[bios_id]] = float(v)
                #print(metabolites[metabolite_id])
            else:
                print('!')
            pass
        #print(object_stoichiometry)
        return object_stoichiometry
    
    def convert_modelreaction(self, r):
        mr_id = r['id']
        name = r['name'] if 'name' in r else r['id']
        lower_bound, upper_bound = (-10,10)#reaction.get_reaction_constraints()
        id = mr_id
        cobra_reaction = Reaction(id, 
                                  name=name, 
                                  lower_bound=lower_bound, 
                                  upper_bound=upper_bound)

        cobra_reaction.add_metabolites(self.convert_modelreaction_stoichiometry(r))

        #gpr = get_gpr(mr)
        #gpr_string = get_gpr_string(gpr)
        #print(gpr_string)
        #reaction.gene_reaction_rule = gpr_string
        return cobra_reaction
    
    def build(self, model_id='model'):
        self.metabolites = {}
        self.reactions = {}
        self.biomass_reactions = set()
        
        compartments = {}
        for o in self.model_cmps:
            if 'name' in o:
                compartments[o['id']] = o['name']
            else:
                compartments[o['id']] = o['id']
            
        for m in self.model_spis:
            self.metabolites[m['bios_id']] = self.convert_modelcompound(m)
            
        for r in self.model_rxns:
            cobra_reaction = self.convert_modelreaction(r)
            if 'Biomass' in cobra_reaction.name:
                self.biomass_reactions.add(cobra_reaction.id)
            #print(cobra_reaction)
            if not cobra_reaction.id in self.reactions:
                self.reactions[cobra_reaction.id] = cobra_reaction
            else:
                print('r111')
        logger.warning(self.biomass_reactions)
        
        cobra_model = Model(model_id)
        cobra_model.compartments = compartments
        cobra_model.add_metabolites(list(self.metabolites.values()))
        cobra_model.add_reactions(list(self.reactions.values()))
        if len(self.biomass_reactions) > 0:
            cobra_model.objective = list(self.biomass_reactions)[0]
        
        return cobra_model
