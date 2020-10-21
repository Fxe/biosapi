from enum import Enum


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


def map_to_identifiers_org_compound(bios_db):
    m = {
        'ModelSeed': Databases.MODELSEED_COMPOUND.value,
        'BiGG': Databases.BIGG1_METABOLITE.value,
        'BiGGMetabolite': Databases.BIGG_METABOLITE.value,
        'LigandCompound': Databases.KEGG_COMPOUND.value,
        'MetaCyc': Databases.METACYC_COMPOUND.value,
    }

    if bios_db in m:
        return m[bios_db]

    return None


def map_to_identifiers_org_reaction(bios_db):
    m = {
        'ModelSeedReaction': Databases.MODELSEED_REACTION.value,
        'BiGG': Databases.BIGG1_REACTION.value,
        'BiGGReaction': Databases.BIGG_REACTION.value,
        'LigandReaction': Databases.KEGG_REACTION.value,
        'MetaCyc': Databases.METACYC_REACTION.value,
    }

    if bios_db in m:
        return m[bios_db]

    return None
