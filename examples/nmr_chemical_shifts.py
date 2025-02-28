import isicle

# Load example structure
geom = isicle.load("CCCC(=O)O")

# Initial optimization
geom = geom.initial_optimize(embed=True, forcefield="UFF", ff_iter=200)

# Molecular dynamics
md = isicle.md.md(
    geom,
    backend="xtb",
    forcefield="gff",
    optlevel="Normal",
    ewin=3,
    task="conformer",
    solvation="water",
    processes=4,
)
# Recommend saving, joblib or pkl
isicle.save("nmr_conformers.joblib", md)

# Parse molecular dynamics
# Now separate from isicle.md.md()
md_parsed = md.parse()

# Using all conformers from CREST
conformers = md_parsed["geometry"]["conformers"]

# Density functional theory
dft = conformers.apply(
    func=isicle.qm.dft,
    backend="nwchem",
    tasks=["energy", "shielding"],
    functional="b3lyp",
    basis_set="6-31g*",
    ao_basis="cartesian",
    solvent="h2o",
    gas=False,
    atoms=["C", "H"],
    temp=298.15,
    scratch_dir="/tmp",
    processes=8,
)

# Recommended saving, joblib or pkl
isicle.save("nmr_dft_conformers.joblib", dft)

# Parse density functional theory results
dft_parsed = [x.parse() for x in dft]

# Create conformational ensemble, which enable func mapping
# Combine shielding result across conformers with Boltzmann weighting
shielding = isicle.conformers.ConformationalEnsemble(
    [x["geometry"] for x in dft_parsed]
).reduce("_shielding")

shifts = isicle.conformers.transform(
    shielding["mean"],
    m={"H": -0.9921, "C": -0.9816},
    b={"H": 32.2773, "C": 180.4295},
    atom=shielding["atom"],
    index=shielding["index"],
)
