import json
import copy
import cobra
from biosapi.bios_model_mapper import BiosModelMapper
from biosapi.io.bios_model_builder import BiosModelToCobraBuilder


def add_suffix(model, suffix):
    model_id = suffix
    model = copy.deepcopy(model)
    for o in model['metabolites']:
        o['id'] += '@' + model_id
    for o in model['reactions']:
        o['id'] += '@' + model_id
        o['metabolites'] = dict(map(
            lambda x: ("{}@{}".format(x[0], model_id), x[1]),
            o['metabolites'].items()))
    return model


def get_cmp_token(cmps):
    if len(cmps) == 1:
        return list(cmps)[0]
    if len(cmps) == 2:
        if 'b' in cmps and 'e' in cmps:
            return 'b'
        if 'e' in cmps and 'c' in cmps:
            return 'c'
        if 'c' in cmps:
            return list(filter(lambda x: not x == 'c', cmps))[0]
    print('!@??!', cmps)
    return None


class MergeModelBuilder:

    def __init__(self, bios):
        self.bios = bios
        self.model_mappers = {}
        self.cobra_models = {}

    def with_model(self, model_id):
        self.model_mappers[model_id] = BiosModelMapper(self.bios, model_id)
        model = BiosModelToCobraBuilder.from_api(model_id, self.bios).build()
        self.cobra_models[model_id] = json.loads(cobra.io.to_json(model))
        return self

    def get_cmp_mapping(self):
        cmp_mapping = {}
        for model_id in self.model_mappers:
            cmp_mapping[model_id] = {}
            mm = self.model_mappers[model_id]
            for o in mm.cmp:
                scmp = o['bios_scmp_entry'] if 'bios_scmp_entry' in o else None
                if scmp is not None:
                    cmp_mapping[model_id][o['id']] = scmp
        return cmp_mapping

    def build(self):
        metabolites = {}
        reactions = {}
        genes = []
        compartments = {}
        cmp_mapping = self.get_cmp_mapping()
        for i in self.model_mappers:
            model = self.cobra_models[i]
            model = add_suffix(model, i)
            mm = self.model_mappers[i]
            spi_mapping = mm.get_spi_annotation('ModelSeed', 5)
            rxn_mapping = mm.get_rxn_annotation('ModelSeedReaction', 4)
            replaced_cpd_ids = {}
            for o in model['metabolites']:
                if 'annotation' not in o:
                    o['annotation'] = {}
                if i not in o['annotation']:
                    o['annotation'][i] = []
                sid = o['id'].split('@')[0]
                if sid not in o['annotation'][i]:
                    o['annotation'][i].append(sid)
                base_id = o['id']
                replaced_cpd_ids[o['id']] = o['id']
                cmp_token = o['compartment']
                if i in cmp_mapping and cmp_token in cmp_mapping[i]:
                    cmp_token = cmp_mapping[i][cmp_token]
                o['compartment'] = cmp_token
                if o['id'].split('@')[0] in spi_mapping:
                    replace_id = spi_mapping[o['id'].split('@')[0]] + '_' + cmp_token
                    replaced_cpd_ids[o['id']] = replace_id
                    o['id'] = replace_id
                else:
                    o['name'] = "{} [{}]".format(o['name'], cmp_token)
                if o['id'] not in metabolites:
                    metabolites[o['id']] = o
                else:
                    if i not in metabolites[o['id']]['annotation']:
                        metabolites[o['id']]['annotation'][i] = []
                    if sid not in metabolites[o['id']]['annotation'][i]:
                        metabolites[o['id']]['annotation'][i].append(sid)
                #metabolites[o['id']]['name'] += ';' + base_id

            for o in model['reactions']:
                base_id = o['id']
                o['metabolites'] = dict(map(
                    lambda x: (replaced_cpd_ids[x[0]], x[1]),
                    o['metabolites'].items()))
                if o['id'].split('@')[0] in rxn_mapping:
                    # print(o['id'], rxn_mapping[o['id'].split('@')[0]])
                    metabolites_cmp = set(map(lambda x: metabolites[x]['compartment'], o['metabolites']))
                    cmp_token = get_cmp_token(metabolites_cmp)
                    # print(o['metabolites'])
                    if cmp_token is not None:
                        o['id'] = rxn_mapping[o['id'].split('@')[0]] + '_' + cmp_token
                    else:
                        print(o['id'], metabolites_cmp)
                if o['id'] not in reactions:
                    reactions[o['id']] = o
                reactions[o['id']]['name'] += ';' + base_id
        merge_model = {
            'metabolites': list(metabolites.values()),
            'reactions': list(reactions.values()),
            'genes': genes,
            'id': 'model',
            'compartments': compartments,
            'version': 'x'
        }
        return merge_model
