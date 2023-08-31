from setuptools import setup, find_packages

setup(
    name="afpairs",
    version="0.1",
    packages=find_packages(),
    author="Christian Poitras",
    author_email="christian.poitras@ircm.qc.ca",
    description="Find interactions using AlphaFold",
    keywords="bioinformatics, AlphaFold",
    url="https://github.com/benoitcoulombelab/alphafold-pairs.git",
    license="GNU General Public License version 3",
    classifiers=[
        "License :: OSI Approved :: GNU General Public License version 3"
    ],
    install_requires=[
        "biopython>=1.81"
    ],
    entry_points={
        "console_scripts": [
            "delete-fasta = afpairs.DeleteFasta:main",
            "fasta-id = afpairs.FastaId:main",
            "fasta-pairs = afpairs.FastaPairs:main",
            "id-convert = afpairs.IdConvert:main",
            "interaction-score = afpairs.InteractionScore:main",
            "multi-interaction-score = afpairs.MultiInteractionScore:main",
            "random-sequences = afpairs.RandomSequences:main",
            "score-matrix = afpairs.ScoreMatrix:main",
            "split-fasta = afpairs.SplitFasta:main",
        ]
    }
)
