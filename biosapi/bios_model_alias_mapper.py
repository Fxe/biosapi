class BiosModelAliasMapper():
    
    def __init__(self, api):
        self.api = api
        pass
    
    def load_model(self, model_id):
        self.model_reactions = self.api.get_model_reactions(model_id)
        self.model_species = self.api.get_model_species(model_id)
        self.model_species_annotation = self.api.get_model_species_annotation(model_id)
        self.model_reactions_annotation = self.api.get_model_reactions_annotation(model_id)
        
        
        self.id_sid = {}
        self.sid_to_spi_id = {}
        self.sid_to_rxn_id = {}
        for o in self.model_species:
            self.id_sid[o['bios_id']] = o['id']
            self.sid_to_spi_id[o['id']] = o['bios_id']
        for o in self.model_reactions:
            self.id_sid[o['bios_id']] = o['id']
            self.sid_to_rxn_id[o['id']] = o['bios_id']
            
    def get_cpd_mapping(self, database, min_score = 0):
        cpd_mapping = {}
        for spi_uid in self.model_species_annotation:
            spi_id = self.id_sid[int(spi_uid)]
            cpd_mapping[spi_id] = set()
            ranked_mapping = {}
            best_score = 0
            for uid in self.model_species_annotation[spi_uid]['annotation']:
                annotation = self.model_species_annotation[spi_uid]['annotation'][uid]
                if annotation['database'] == database:
                    score = annotation['link']['authors']
                    score = dict(map(lambda x : x.split(':'), score.split(';')))
                    max_score = 0
                    for usr in score:
                        if int(score[usr]) > max_score:
                            max_score = int(score[usr])
                    if max_score > best_score:
                        best_score = max_score
                    if not max_score in ranked_mapping:
                        ranked_mapping[max_score] = set()
                    ranked_mapping[max_score].add(annotation['entry'])
                    
            if len(ranked_mapping) > 0 and best_score >= min_score:
                for i in ranked_mapping[best_score]:
                    cpd_mapping[spi_id].add(i)
        return cpd_mapping
    
    def get_rxn_mapping(self, database, min_score = 0):
        rxn_mapping = {}
        for rxn_uid in self.model_reactions_annotation:
            rxn_id = self.id_sid[int(rxn_uid)]
            rxn_mapping[rxn_id] = set()
            for uid in self.model_reactions_annotation[rxn_uid]['annotation']:
                annotation = self.model_reactions_annotation[rxn_uid]['annotation'][uid]
                if annotation['database'] == database:
                    rxn_mapping[rxn_id].add(annotation['entry'])
        return rxn_mapping