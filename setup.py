from setuptools import setup, find_packages

with open("isicle/__init__.py") as f:
    exec([x for x in f.readlines() if "__version__" in x][0])

with open("README.md") as f:
    readme = f.read()

with open("LICENSE") as f:
    license = f.read()

with open("requirements.txt", "r") as f:
    install_requires = f.read().splitlines()

pkgs = find_packages(exclude=("examples", "docs", "resources"))

setup(
    name="isicle",
    version=__version__,  # TODO: switch to versioneer.get_version() or alternate method
    description="ISiCLE: in silico chemical library engine",
    long_description=readme,
    long_description_content_type="text/markdown",
    author="Sean M. Colby",
    author_email="sean.colby@pnnl.gov",
    license=license,
    packages=pkgs,
    python_requires=">=3.9",
    install_requires=install_requires,
    classifiers=[
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
)
