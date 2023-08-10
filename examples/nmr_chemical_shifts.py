import isicle

# Load example structure
geom = isicle.load('CCCC(=O)O')

# Initial optimization
geom = geom.initial_optimize(embed=True, forcefield="UFF", ff_iter=200)

# Molecular dynamics
md_result = geom.md(program='xtb',
                    task='conformer',
                    forcefield='gff',
                    ewin=1,
                    optlevel='Normal')

# Density functional theory
dft_result = md_result.get_structures().apply(func=isicle.qm.dft,
                                              tasks=['energy', 'shielding'],
                                              functional='b3lyp',
                                              basis_set='3-21g*',
                                              ao_basis='cartesian',
                                              charge=0,
                                              atoms=['C', 'H'],
                                              temp=298.15,
                                              processes=8)

# Combine shielding result across conformers
shielding = dft_result.get_structures().reduce('shielding', func='boltzmann')
