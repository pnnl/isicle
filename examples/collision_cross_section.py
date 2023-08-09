import isicle

# Load example structure
geom = isicle.load('CCCC(=O)O')

# Perform initial optimization
geom = geom.initial_optimize(embed=True, forcefield="UFF", ff_iter=200)

# Form adducts
ionization_res = geom.ionize(method='explicit',
                             ion_list=['H-'],
                             element_list=['O'])

# Iterate through adduct ions
ccs_container = []
for adduct in ionization_res.get_structures():
    # Molecular dynamics
    md_result = adduct.md(program='xtb',
                          task='conformer',
                          forcefield='gff',
                          ewin=1,
                          optlevel='Normal',
                          charge=adduct.charge)

    # Density functional theory
    dft_result = md_result.get_structures().apply(func=isicle.qm.dft,
                                                  tasks=['optimize', 'energy'],
                                                  functional='b3lyp',
                                                  basis_set='3-21g*',
                                                  ao_basis='cartesian',
                                                  charge=adduct.charge,
                                                  temp=298.15,
                                                  processes=8)

    # Calculate CCS
    ccs_result = dft_result.get_structures().apply(isicle.mobility.calculate_ccs,
                                                   lennard_jones='default',
                                                   i2=5013489,
                                                   buffer_gas='nitrogen',
                                                   buffer_gas_mass=28.014,
                                                   temp=300,
                                                   ipr=1000,
                                                   itn=10,
                                                   inp=48,
                                                   imp=1024,
                                                   processes=8)

    # Combine shielding result across conformers
    ccs = ccs_result.get_structures().reduce('ccs', func='boltzmann')

    # Add to container
    ccs_container.append(ccs)
