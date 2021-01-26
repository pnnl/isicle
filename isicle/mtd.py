
def xtb(input, method = 'gff', optlevel='Normal'):
	"""
	Geometry input can be TURBOMOLE coordinates (.coord/.tmol), any valid Xmol (e.g. .xyz), 
	mol files (.mol), Structure-Data files (.sdf), Protein Database Files (.pdb), 
	Vasp’s POSCAR and CONTCAR files (.poscar/.contcar/.vasp) and DFTB+ genFormat files (.gen)
	"""
	# TODO: define output
	raise NotImplementedError

def crest():
    # TODO: define input, arguments, and output
    # Should return instance (or list) of CREST conformers
    raise NotImplementedError

def enso():
	# TODO: define input, arguments, and output
	# Should return instance (or list) of ENSO calculated chemical shifts
	raise NotImplementedError

def amber():
    # Don't work on this one yet
    # TODO: define input, arguments, and output
    # Should return instance (or list) of MDOptimizedGeometry
    raise NotImplementedError

# TODO: add other MD methods as needed
