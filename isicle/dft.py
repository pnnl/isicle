"""
Assumes NWChem as default quantum chemistry package. Other options can be added later:
 Orca, Quantum Espresso, Gaussian, Schrodinger, etc.
"""



def optimize(self, program= 'NWChem', functional='b3lyp', basisset='6-31g*'):
	"""
	Optimize geometry, either XYZ or PDB, using stated functional and basis set.
	Additional inputs can be grid size, optimization criteria level, 
	"""
    # TODO: define input, arguments, and output
    # Should return instance (or list) of DFT
    raise NotImplementedError

def property():
	"""
	Calculates a property (shielding tensor, frequency, charges) of a given geometry
	"""
	# TO DOL define input, arguments, output
	# Should return list of freq/shifts/charges of a given property
	raise NotImplementedError

# TODO: add other DFT methods as needed


