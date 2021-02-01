from isicle.interfaces import AdductInterface
from isicle.geometry import Geometry


class Adduct(Geometry, AdductInterface):

    def __init__(self):
        self.geom = None
        self.ions = None
        self.anions = None
        self.cations = None
        self.complex = None

    def load_geometry(self, path: str):
        '''
        Reads geometry3D object from file.
        '''
        self.geom = load(path)
        return self.geom

    def load_ions(self, path: str):
        '''
        Reads adduct ions from file.
        '''
        # load text file
        self.ions = load(path).split(',')
        return self.ions

    def parse_ions(self):
        '''
        Changes order of ionization from specification in ion list file.
        '''
        anions = []
        cations = []
        complex = []

        for x in ions:
            if ('+' and '-') in x:
                complex.append(x)
            elif ('+') in x:
                cations.append(x)
            elif ('-') in x:
                anions.append(x)
            else:
                raise Exception('Unrecognized ion specified.')
                # TODO insert more descriptive error

        self.anions = anions
        self.cations = cations
        self.complex = complex
        return self.anions, self.cations, self.complex

    def generate_adducts(self, alt=False, inplace=False):
        '''
        Calls different generation methods according to desired ion reaction.
        '''
        # TODO check if desired adduct can be formed
        self.check_valid()
        # TODO institute iterative approach to complex adduct formation
        # consult https://github.com/grimme-lab/crest/issues/2 for supported ions by xtb CREST
        if alt == True:
            self.crest_generation()
        else:
            self.manual_generation()

        # TODO institute manual generation capabilities of output molecule
        # default positive_mode/negative_mode to xtb CREST, alt for manual generation

    def check_valid(self, ion):
        '''
        Performs substructure search and checks if proposed adduct can be formed.
        '''
        # consider calc_new_mz from
        # https://github.com/pnnl/mame/blob/24f9fcb19639d5a1c2ca0f788c55d6a2efe18ca6/mame/mame_utils.py
        # consider add_formala, remove_formula from
        # https://github.com/pnnl/mame/blob/24f9fcb19639d5a1c2ca0f788c55d6a2efe18ca6/mame/formula_module.py
        mol.HasSubstructMatch((ion.strip('+')).strip('-'))
        # complex structure checking, spectre
        # consider how substructure would break off, not all atoms may be from same region

    def save(self):
        '''
        Writes output molecules to file.
        '''
        # pull from geometry

    def manual_genertion(self):
        '''
        Creates specified adducts at every atom site.
        '''
        # TODO iterative method for adduct generation
        # for adduct type specified in parse ions, call positive mode or negative mode
        ion_lst = self.ions

    def manual_positive_mode(self):
        '''
        Adds cation to 3D molecular structure.
        '''
        self.add_ion()

    def manual_negative_mode(self):
        '''
        Pulls anion from 3D molecular structure.
        '''
        self.remove_ion()

    def _remove_ion():
        '''
        '''

    def remove_ion():
        '''
        '''
        self._remove_ion()

    def _add_ion():
        '''
        '''

    def add_ion():
        '''
        '''
        self._add_ion()

    def crest_generation(self):
        '''
        Alternative adduct generation method.
        Utilizes xtb CREST software from Grimme group.
        Documentation:
        '''
        self.crest_positive_mode()
        self.crest_negative_mode()

    def crest_positive_mode(self, path, ion=None, alt=False):
        '''
        Runs xtb CREST protonate tool
        https://xtb-docs.readthedocs.io/en/latest/crestxmpl.html#protonation-site-screening

        Default:
        [H+]: (ex.) crest 'struc.xyz' -protonate

        Alternate:
        (ex.) crest 'struc.xyz' -protonate -swel STR
        STR must be string of element and charge.
        [Na+]: (ex.) crest 'struc.xyz' -protonate -swel 'Na+'
        '''
        # consider -twel for tautomeric forms of output adduct
        if alt:
            # crest()
            # "crest %s -protonate -swel %s"%(path,ion)
        else:
            # crest()
            # "crest %s -protonate"%(path)

    def crest_negative_mode(self, path, ion=None, alt=False):
        '''
        Runs xtb CREST deprotonate tool

        https://xtb-docs.readthedocs.io/en/latest/crestxmpl.html#deprotonation-site-screening

        Default:
        [H-]: (ex.) crest 'struc.xyz' -deprotonate
        '''
        # consider -twel for tautomeric forms of output adduct
        if alt:
            # crest()
            # "crest %s -deprotonate -swel %s"%(path,ion)
        else:
            # crest()
            # "crest %s -deprotonate"%(path)
