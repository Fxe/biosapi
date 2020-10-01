

class BiosModelMapper:
    
    def __init__(self, bios, model_id):
        self.model_id = model_id
        self.cmp = bios.get_model_compartments(model_id)
        self.rxn = {}
        self.spi = {}
        for o in bios.get_model_species(model_id):
            self.spi[o['bios_id']] = o
        for o in bios.get_model_reactions(model_id):
            self.rxn[o['bios_id']] = o
        self.rxn_annotation = bios.get_model_reactions_annotation(model_id)
        self.spi_annotation = bios.get_model_species_annotation(model_id)
        self.spi_sid_to_uid = {}
        self.rxn_sid_to_uid = {}
        
        for uid in self.spi:
            o = self.spi[uid]
            if 'id' in o and not o['id'] in self.spi_sid_to_uid:
                self.spi_sid_to_uid[o['id']] = uid
            else:
                print('!')
        for uid in self.rxn:
            o = self.rxn[uid]
            if 'id' in o and not o['id'] in self.rxn_sid_to_uid:
                self.rxn_sid_to_uid[o['id']] = uid
            else:
                print('!')

    def get_rxn_annotation(self, database, min_score):
        result = {}
        for mrxn_uid in self.rxn_annotation:
            annotation = self.rxn_annotation[mrxn_uid]['annotation']
            for rxn_uid in annotation:
                rxn_annotation = annotation[rxn_uid]
                if rxn_annotation['database'] == database:
                    score = dict(map(lambda x: x.split(':'), rxn_annotation['link']['authors'].split(';')))
                    for user_id in score:
                        if int(score[user_id]) >= min_score:
                            result[self.rxn[int(mrxn_uid)]['id']] = rxn_annotation['entry']

        return result

    def get_spi_annotation(self, database, min_score):
        result = {}
        for spi_uid in self.spi_annotation:
            annotation = self.spi_annotation[spi_uid]['annotation']
            for cpd_uid in annotation:
                cpd_annotation = annotation[cpd_uid]
                if cpd_annotation['database'] == database:
                    score = dict(map(lambda x: x.split(':'), cpd_annotation['link']['authors'].split(';')))
                    for user_id in score:
                        if int(score[user_id]) >= min_score:
                            result[self.spi[int(spi_uid)]['id']] = cpd_annotation['entry']

        return result

    def get_spi_database_mapping(self, cpd_id, database):
        match_uid = set()
        for spi_uid in self.spi_annotation:
            o = self.spi_annotation[spi_uid]['annotation']
            for cpd_uid in o:
                if o[cpd_uid]['database'] == database and o[cpd_uid]['entry'] == cpd_id:
                    match_uid.add(int(spi_uid))

        match_spi = set()
        for spi_uid in match_uid:
            match_spi.add(self.spi[spi_uid]['id'])
        return match_spi
    
    def get_rxn_database_mapping(self, rxn_id, database):
        match_uid = set()
        for mrxn_uid in self.rxn_annotation:
            o = self.rxn_annotation[mrxn_uid]['annotation']
            for rxn_uid in o:
                if o[rxn_uid]['database'] == database and o[rxn_uid]['entry'] == rxn_id:
                    match_uid.add(int(mrxn_uid))

        match_rxn = set()
        for mrxn_uid in match_uid:
            match_rxn.add(self.rxn[mrxn_uid]['id'])
        return match_rxn
                
    def get_spi_uid(self, sid):
        if sid in self.spi_sid_to_uid:
            return self.spi_sid_to_uid[sid]
        return None
    
    def get_rxn_uid(self, sid):
        if sid in self.rxn_sid_to_uid:
            return self.rxn_sid_to_uid[sid]
        return None
