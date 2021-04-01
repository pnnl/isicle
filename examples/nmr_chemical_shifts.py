import isicle

# Load example structure
geom = isicle.geometry.load('example.smi')

# Optional desalt
geom = geom.desalt()

# Optional neutralize
geom = geom.neutralize()

# Optional tautomerize
geom = geom.tautomerize()

# Molecular dynamics
kwargs = {}
conformers, md_result = geom.md_optimize(program='xtb',
                                         **kwargs)

# Density functional theory
dft_result = conformers.apply(func=isicle.qm.dft,
                              tasks=['optimize', 'shielding'],
                              functional='b3lyp',
                              basis_set='6-31g*',
                              ao_basis='cartesian',
                              charge=0,
                              atoms=['C', 'H'],
                              frequency=True,
                              temp=298.15)

# Separate conformers from result
conformers = isicle.conformers.build_conformational_ensemble(
    [x[0] for x in dft_result])
dft_result = [x[1] for x in dft_result]

# Combine shielding result across conformers
shielding = conformers.reduce('shielding', func='boltzmann')

# Convert to shifts
# TODO: implement shielding to shifts conversion
