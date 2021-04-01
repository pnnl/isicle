import isicle

# Load example structure
geom = isicle.geometry.load('example.smi')

# Desalt
geom = geom.desalt()

# Neutralize
geom = geom.neutralize()

# Tautomerize
geom = geom.tautomerize()

# Form adducts
adducts = geom.generate_adducts()

ccs_result = {}
for label, adduct in adducts.items():
    # Molecular dynamics
    kwargs = {}
    conformers, md_result = geom.md_optimize(program='xtb',
                                             **kwargs)

    # Density functional theory
    dft_result = conformers.apply(func=isicle.qm.dft,
                                  tasks=['optimize'],
                                  functional='b3lyp',
                                  basis_set='6-31g*',
                                  ao_basis='cartesian',
                                  charge=adduct.charge,  # Need to resolve this
                                  frequency=True,
                                  temp=298.15)

    # Separate conformers from result
    conformers_opt = isicle.conformers.build_conformational_ensemble(
        [x[0] for x in dft_result])
    dft_result = [x[1] for x in dft_result]

    # Combine shielding result across conformers
    ccs = conformers_opt.reduce('shielding', func='boltzmann')

    # Add to dictionary
    ccs_result[label] = ccs
