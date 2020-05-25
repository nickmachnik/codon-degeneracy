import os
from setuptools import setup, find_packages, Command


__version__ = None
exec(open('version.py').read())


class CleanCommand(Command):
    """
    Custom clean command to tidy up the project root.
    """
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        os.system(
            'rm -vrf ./build ./dist ./*.pyc ./*.tgz ./*.egg-info ./htmlcov')


setup(
    name='codon-degeneracy',
    version=__version__,
    description="""
        Routines for the extraction of degenerate sides and calculation of rates
        of neutral substitutions from sequences and alignments.""",
    # long_description_content_type="text/markdown",
    url='https://github.com/nickmachnik/codon-degeneracy',
    setup_requires=[
        'setuptools>=18.0',
    ],
    packages=find_packages(),
    install_requires=[
        'scikit-bio',
        'biopython'
    ],
    # scripts=['bin/chess'],
    cmdclass={
        'clean': CleanCommand
    },
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
)