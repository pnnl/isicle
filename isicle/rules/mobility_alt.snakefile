from os.path import *
from isicle.utils import getOS
from pkg_resources import resource_filename


# snakemake configuration
include: 'adducts.snakefile'
IMPACT = resource_filename('isicle', 'resources/IMPACT/%s/impact' % getOS())


rule impact:
    input:
        rules.touchAdducts.output.xyz
    output:
        abspath(join('output', 'mobility', 'impact', 'runs', '{id}_{adduct}_{addID}.txt'))
    version:
        '{IMPACT} -version -nocite'
    log:
        abspath(join('output', 'mobility', 'impact', 'runs', 'logs', '{id}_{adduct}_{addID}.log'))
    benchmark:
        abspath(join('output', 'mobility', 'impact', 'runs', 'benchmarks', '{id}_{adduct}_{addID}.benchmark'))
    # group:
    #     'mobility_alt'
    shell:
        # run impact on adducts
        'IMPACT_RANDSEED={config[impact][seed]} {IMPACT} {input} -o {output} -H \
         -shotsPerRot {config[impact][shotsPerRot]} -convergence {config[impact][convergence]} \
         -nRuns {config[impact][nRuns]} -nocite &> {log}'


rule postprocess:
    input:
        ccs = rules.impact.output,
        mass = rules.calculateMass.output
    output:
        he = abspath(join('output', 'mobility', 'impact', 'ccs', '{id}_{adduct}_{addID}.He.ccs')),
        n2 = abspath(join('output', 'mobility', 'impact', 'ccs', '{id}_{adduct}_{addID}.N2.ccs'))
    version:
        'isicle --version'
    log:
        abspath(join('output', 'mobility', 'impact', 'ccs', 'logs', '{id}_{adduct}_{addID}.benchmark'))
    benchmark:
        abspath(join('output', 'mobility', 'impact', 'ccs', 'benchmarks', '{id}_{adduct}_{addID}.benchmark'))
    # group:
    #     'mobility_alt'
    shell:
        'python -m isicle.scripts.parse_impact {input.ccs} {input.mass} {output.he} {output.n2} \
         --alpha {config[ccs][alpha]} --beta {config[ccs][beta]} &> {log}'
