#!/usr/bin/env nextflow

params.cwd = System.getProperty("user.dir")

process conformers {
  label 'conformers'
  input:
    path inputfile
  output:
    path '*.joblib'
    publishDir 'output/conformers', mode: 'copy'
  """
  #!/usr/bin/env python3
  import isicle
  import yaml
  import os
  from os.path import *

  cwd = "${params.cwd}"

  with open(f"{cwd}/config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

  geom = isicle.io.load("${inputfile}")
  geom = geom.initial_optimize(embed=True)
  charge = geom.get_formal_charge()

  conformers = geom.md(forcefield=config['md']['forcefield'],
                      ewin=config['md']['energywindow'],
                      task=config['md']['task'],
                      charge=charge,
                      processes=config['md']['threads'])

  for conformer in conformers.geom:
    conformer.set_formal_charge(charge)
    conformer.basename = geom.basename
    output_path = join(f'{geom.basename}_{conformer.conformerID}_md.joblib')
    isicle.io.save(output_path, conformer)
  """
}

process dft {
  label 'dft'
  input:
    path inputfile
  output:
    path '*.joblib'
    publishDir 'output/dft', mode: 'copy'
  """
  #!/usr/bin/env python3
  import isicle
  import yaml
  import os
  from os.path import *

  cwd = "${params.cwd}"

  with open(f"{cwd}/config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

  geom = isicle.io.load("${inputfile}")

  dft = isicle.qm.dft(geom, tasks=config['dft']['tasks'],
                      functional=config['dft']['functional'],
                      basis_set=config['dft']['basis_set'],
                      ao_basis=config['dft']['ao_basis'],
                      charge=geom.charge,
                      atoms=config['dft']['atoms'],
                      temp=config['dft']['temp'],
                      cosmo=config['dft']['cosmo'],
                      solvent=config['dft']['solvent'],
                      gas=config['dft']['gas'],
                      max_iter=config['dft']['max_iter'],
                      mem_global=config['dft']['mem_global'],
                      mem_heap=config['dft']['mem_heap'],
                      mem_stack=config['dft']['mem_stack'],
                      scratch_dir=config['dft']['scratch_dir'],
                      processes=4,
                      command=config['dft']['command'])

  output_path = join(f'{geom.basename}_{geom.conformerID}_dft.joblib')
  dft.basename = geom.basename
  dft.conformerID = geom.conformerID
  isicle.io.save(output_path, dft)
  """
}

workflow {
  data = channel.fromPath('input/*.*') | flatten | conformers | flatten | dft 
}
