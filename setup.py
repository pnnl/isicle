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
    package_data={'': ['rules/*.snakefile',
                       'resources/amber/leaprc.ff14SB',
                       'resources/amber/leaprc.gaff',
                       'resources/amber/*.template',
                       'resources/mobcal/mobcal_cascade',
                       'resources/mobcal/mobcal_constance',
                       'resources/mobcal/*.params',
                       'resources/mobcal/atomic_mass.tsv',
                       'resources/mobcal/atomtype_parameters.in',
                       'resources/IMPACT/linux/impact',
                       'resources/IMPACT/osx/impact',
                       'resources/nwchem/*.smi',
                       'resources/nwchem/*.template']},
    include_package_data=True
)
