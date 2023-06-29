#!/usr/bin/env nextflow

params.cwd = System.getProperty("user.dir") 

process adducts {
  input:
    file inputfile
  output:
    path '*.joblib'
    publishDir 'output/adducts/', mode: 'copy'
  """
  #!/usr/bin/env python3
  import isicle
  import yaml
  import os
  import sys
  from os.path import *

  
  cwd = "${params.cwd}"
  with open(f"{cwd}/config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

  geom = isicle.io.load("${inputfile}")
  geom = geom.initial_optimize(embed=True)

  adducts = geom.ionize(ion_list=config['adducts']['ion_list'],
                        element_list=config['adducts']['element_list'],
                        method=config['adducts']['method'])
            
  for add in adducts.adducts:
    output_path = f'{add.basename}_{add.ion}_{add.adductID}.joblib'
    isicle.io.save(output_path, add)
  """
}

process conformers {
  input:
    file inputfile
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

  add = isicle.io.load("${inputfile}")
  charge = add.get_formal_charge()

  conformers = add.md(forcefield=config['md']['forcefield'],
                      ewin=config['md']['energywindow'],
                      task=config['md']['task'],
                      charge=charge,
                      processes=config['md']['threads'])


  for conformer in conformers.geom:
    conformer.set_formal_charge(charge)
    conformer.basename = add.basename
    conformer.ion = add.ion
    conformer.adductID = add.adductID
    output_path = join(f'{add.basename}_{add.ion}_{add.adductID}_{conformer.conformerID}_md.joblib')
    isicle.io.save(output_path, conformer)
  """
}

process dft {
  input:
    file inputfile
  output:
    path '*.joblib'
    publishdir 'output/dft', mode: 'copy'
  """
  #!/usr/bin/env python3
  import isicle
  import yaml
  import os
  from os.path import *

  cwd = "${params.cwd}"

  with open(f"{cwd}/config.yaml", "r") as f:
    config = yaml.load(f, loader=yaml.fullloader)

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
                      processes=config['dft']['processes'])

  dft.basename = geom.basename
  dft.conformerid = geom.conformerid
  dft.adductid = geom.adductid
  output_path = f'{geom.basename}_{geom.ion}_{geom.adductid}_{geom.conformerid}_dft.joblib'
  isicle.io.save(output_path, dft)
  """
}

// Onlt available on linux
process mobility {
  input:
    file inputfile
  output:
    path '*.joblib'
    publishdir 'output/mobility', mode: 'copy'
  """
  #!/usr/bin/env python3
  import isicle
  import yaml
  import os
  from os.path import *

  cwd = "${params.cwd}"

  with open(f"{cwd}/config.yaml", "r") as f:
    config = yaml.load(f, loader=yaml.fullloader)

  geom = isicle.io.load("${inputfile}")

  ccs = isicle.mobility.calculate_ccs(geom.geom, lennard_jones=config['mobility']['lennard_jones'],
                                      i2=config['mobility']['i2'],
                                      buffer_gas=config['mobility']['buffer_gas'],
                                      buffer_gas_mass=config['mobility']['buffer_gas_mass'],
                                      temp=config['mobility']['temp'],
                                      ipr=config['mobility']['ipr'],
                                      itn=config['mobility']['itn'],
                                      inp=config['mobility']['inp'],
                                      imp=config['mobility']['imp'],
                                      command=config['mobility']['command'])

  ccs.basename = geom.basename
  ccs.conformerid = geom.conformerid
  ccs.adductid = geom.adductid
  ccs.ion = geom.ion
  output_path = f'{geom.basename}_{geom.ion}_{geom.adductid}_{geom.conformerid}_dft.joblib'
  isicle.io.save(output_path, ccs)
  """
}
workflow {
  data = channel.fromPath('input/*') | flatten | adducts | flatten | conformers | flatten | dft | mobility}
