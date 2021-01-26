from isicle.interfaces import AdductInterface


class Adduct(AdductInterface):

    def __init__(self):
        self.contents = None

    def load(self):
        """
        Loads geometry3D object.
        """

    def get_ions(self):
        """
        Parses in adduct ions to be generated
        """
        # parse in ions from file
        # self.ions =

    def generate_adducts(self):
        """
        Calls different generation methods according to desired ion reaction.
        """
        # TODO check if desired adduct can be formed
        self.check_valid()
        # TODO institute iterative approach to complex adduct formation
        # consult https://github.com/grimme-lab/crest/issues/2 for supported ions by xtb CREST
        self.positive_mode()
        self.negative_mode()
        # TODO institute manual generation capabilities of output molecule
        # default positive_mode/negative_mode to xtb CREST, alt for manual generation

    def check_valid(self):
        """
        Checks if proposed adduct can be formed.
        """
        # consider calc_new_mz from
        # https://github.com/pnnl/mame/blob/24f9fcb19639d5a1c2ca0f788c55d6a2efe18ca6/mame/mame_utils.py
        # consider add_formala, remove_formula from
        # https://github.com/pnnl/mame/blob/24f9fcb19639d5a1c2ca0f788c55d6a2efe18ca6/mame/formula_module.py

    def save(self):
        """
        Writes output molecules to file.
        """

    def positive_mode(self):
        """
        Runs xtb CREST protonate tool

        https://xtb-docs.readthedocs.io/en/latest/crestxmpl.html#protonation-site-screening

        Default:
        [H+]: (ex.) crest 'struc.xyz' -protonate

        Alternate:
        (ex.) crest 'struc.xyz' -protonate -swel STR
        STR must be string of element and charge.
        [Na+]: (ex.) crest 'struc.xyz' -protonate -swel 'Na+'
        """
        # consider -twel for tautomeric forms of output adduct

    def negative_mode(self):
        """
        Runs xtb CREST deprotonate tool

        https://xtb-docs.readthedocs.io/en/latest/crestxmpl.html#deprotonation-site-screening

        Default:
        [H-]: (ex.) crest 'struc.xyz' -deprotonate
        """
        # consider -twel for tautomeric forms of output adduct
        # call
