import os
import sys
import subprocess


try:
    import setuptools
except ImportError:
    sys.exit("setuptools package not found. "
             "Please use 'pip install setuptools' first")

from setuptools import setup
from distutils.command.build import build as DistutilsBuild
from distutils.spawn import find_executable

from sd.__version__ import __version__


# Make sure we're running from the setup.py directory.
script_dir = os.path.dirname(os.path.realpath(__file__))
if script_dir != os.getcwd():
    os.chdir(script_dir)


requirements_fn = os.path.join(script_dir, 'requirements.txt')
requirements = []
with open(requirements_fn) as f:
    for line in f:
        line = line.strip()
        if line == 'python-edlib':
            requirements.append('edlib')
        else:
            requirements.append(line)


description = \
    """
StringDecomposer (SD) algorithm takes the set of monomers
and a long error-prone read (or a genomic segment)
and partitions this read into distinct monomers,
providing an accurate translation of each read
from a nucleotide alphabet into a monomer alphabet.
"""


class MakeBuild(DistutilsBuild):
    def run(self):
        if not find_executable("make"):
            sys.exit("ERROR: 'make' command is unavailable")
        try:
            subprocess.check_call(["make"])
        except subprocess.CalledProcessError as e:
            sys.exit("Compilation error: ", e)
        DistutilsBuild.run(self)


setup(
    name="StringDecomposer",
    version=__version__,
    description=description,
    url='https://github.com/ablab/stringdecomposer',
    author='Tatiana Dvorkina',
    author_email='tanunia@gmail.com',
    license='BSD-3-Clause',
    install_requires=requirements,
    packages=['sd'],
    package_dir={'sd': 'sd'},
    package_data={'sd': ['**/*']},
    entry_points={
        'console_scripts': ['run_decomposer=sd.run_decomposer:main']
    },
    cmdclass={'build': MakeBuild}
)
