from enum import Enum


class SubcellularCompartment(Enum):
    
    ANY = 'z'
    MULTIPLE = '*'
    UNKNOWN = '?'
    
    BOUNDARY = 'b'
    EXTRACELLULAR = 'e'
    CYTOSOL = 'c'
    MITOCHONDRIA = 'm'
    VACUOLE = 'v'
    NUCLEUS = 'n'
    GOLGI = 'g'
    RETICULUM = 'r'
    
    PERIPLASM = 'p'
    PEROXISOME = 'x'
    LIPID = 'l'
  
  #LYSOSOME,
  #GLYOXYSOME,
  #CARBOXYSOME,
 # 
  #THYLAKOID_LUMEN,
  #THYLAKOID_MEMBRANE,
 # 
 #        GOLGI_MEMBRANE,
 # MITOCHONDRIA_MEMBRANE,
 #     NUCLEAR_MEMBRANE,
 #      PLASMA_MEMBRANE, //wth ?
 #    VACUOLAR_MEMBRANE,
 # PEROXISOMAL_MEMBRANE,
 #   RETICULUM_MEMBRANE,
 # CELL_ENVELOPE,
 # 
 # CELL_WALL,
 # PLASTID,


scmp_m = {
    'c': SubcellularCompartment.CYTOSOL,
    'm': SubcellularCompartment.MITOCHONDRIA,
    'e': SubcellularCompartment.EXTRACELLULAR,
    'x': SubcellularCompartment.PEROXISOME,
    'v': SubcellularCompartment.VACUOLE,
    'r': SubcellularCompartment.RETICULUM,
    'n': SubcellularCompartment.NUCLEUS,
    'g': SubcellularCompartment.GOLGI,
    'b': SubcellularCompartment.BOUNDARY,
}
