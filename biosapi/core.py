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