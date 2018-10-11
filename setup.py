from setuptools import setup, find_packages


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

with open('requirements.txt') as f:
    required = f.read().splitlines()
    required = None

pkgs = find_packages(exclude=('examples', 'docs'))

setup(
    name='ISiCLE',
    version='0.1.0',
    description='In Silico Chemical Library Engine.',
    long_description=readme,
    author='Sean M. Colby',
    author_email='sean.colby@pnnl.gov',
    url='https://github.com/pnnl/isicle',
    license=license,
    packages=pkgs,
    install_requires=required,
    entry_points={
        'console_scripts': ['isicle = isicle.isicle:cli']
    }
)
