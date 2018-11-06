from os.path import *
from isicle.utils import getOS
from pkg_resources import resource_filename


# snakemake configuration
include: 'adducts.snakefile'
IMPACT = resource_filename('isicle', 'resources/IMPACT/%s/impact' % getOS())


rule impact:
    input:
        rules.generateAdduct.output.xyz
    output:
        join('output', 'mobility', 'impact', 'runs', '{id}_{adduct}.txt')
    version:
        '{IMPACT} -version -nocite'
    log:
        join('output', 'mobility', 'impact', 'runs', 'logs', '{id}_{adduct}.log')
    benchmark:
        join('output', 'mobility', 'impact', 'runs', 'benchmarks', '{id}_{adduct}.benchmark')
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
        he = join('output', 'mobility', 'impact', 'ccs', '{id}_{adduct}.He.ccs'),
        n2 = join('output', 'mobility', 'impact', 'ccs', '{id}_{adduct}.N2.ccs')
    version:
        'isicle --version'
    log:
        join('output', 'mobility', 'impact', 'ccs', 'logs', '{id}_{adduct}.benchmark')
    benchmark:
        join('output', 'mobility', 'impact', 'ccs', 'benchmarks', '{id}_{adduct}.benchmark')
    # group:
    #     'mobility_alt'
    shell:
        'python -m isicle.scripts.parse_impact {input.ccs} {input.mass} {output.he} {output.n2} \
         --alpha {config[ccs][alpha]} --beta {config[ccs][beta]} &> {log}'
