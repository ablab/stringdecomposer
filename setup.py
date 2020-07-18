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
    install_requires=['biopython', 'edlib', 'joblib', 'numpy', 'pandas',
                      'setuptools'],
    packages=['sd'],
    package_dir={'sd': 'sd'},
    package_data={'sd': ['**/*']},
    entry_points={
        'console_scripts': ['sd=sd.run_decomposer:main',
                            'sd_cluster_sequences=sd.cluster_sequences:main']
    },
    cmdclass={'build': MakeBuild}
)
