'''
Assumes NWChem as default quantum chemistry package. Other options can be added later:
 Orca, Quantum Espresso, Gaussian, Schrodinger, etc.
'''


def dft(self, program='NWChem', functional='b3lyp', basisset='6-31g*'):
    '''
    Optimize geometry, either XYZ or PDB, using stated functional and basis set.
    Additional inputs can be grid size, optimization criteria level,
    '''
    # TODO: define input, arguments, and output
    # Should return instance (or list) of DFT

    # save geometry to input file
    # load/generate .nw script
    # submit job
    # parse result
    # return decorated `*OptimizedGeometry` class instance
    raise NotImplementedError

# TODO: add other DFT methods as needed
