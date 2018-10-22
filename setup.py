from setuptools import setup, find_packages
from isicle import __version__


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

with open('requirements.txt') as f:
    required = f.read().splitlines()
    required = None

pkgs = find_packages(exclude=('examples', 'docs', 'resources'))

setup(
    name='isicle',
    version=__version__,
    description='ISiCLE: in silico chemical library engine',
    long_description=readme,
    author='Sean M. Colby',
    author_email='sean.colby@pnnl.gov',
    url='https://github.com/pnnl/isicle',
    license=license,
    packages=pkgs,
    install_requires=required,
    entry_points={
        'console_scripts': ['isicle = isicle.isicle:cli']
    },
    package_data={'isicle': ['isicle/rules/*.snakefile',
                             'isicle/resources/amber/leaprc.ff14SB',
                             'isicle/resources/amber/leaprc.gaff',
                             'isicle/resources/amber/*.template',
                             'isicle/resources/mobcal/mobcal_cascade',
                             'isicle/resources/mobcal/mobcal_constance',
                             'isicle/resources/mobcal/*.params',
                             'isicle/resources/mobcal/atomic_mass.tsv',
                             'isicle/resources/mobcal/atomtype_parameters.in',
                             'isicle/resources/IMPACT/linux/impact',
                             'isicle/resources/IMPACT/osx/impact',
                             'isicle/resources/nwchem/*.smi']},
    include_package_data=True
)
