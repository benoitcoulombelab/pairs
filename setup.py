from setuptools import setup, find_packages

setup(
    name="pairs",
    version="0.1",
    packages=find_packages(),
    author="Christian Poitras",
    author_email="christian.poitras@ircm.qc.ca",
    description="Find interactions using AlphaFold",
    keywords="bioinformatics, AlphaFold",
    url="https://github.com/benoitcoulombelab/pairs.git",
    license="GNU General Public License version 3",
    classifiers=[
        "License :: OSI Approved :: GNU General Public License version 3"
    ],
    install_requires=[
        "biopython>=1.81",
        "pandas>=2.1.0",
        "smokesignal>=0.7",
        "tqdm>=4.66.1"
    ],
    entry_points={
        "console_scripts": [
            "consensus-interface = pairs.ConsensusInterface:main",
            "delete-fasta = pairs.DeleteFasta:main",
            "fasta-id = pairs.FastaId:main",
            "fasta-pairs = pairs.FastaPairs:main",
            "id-convert = pairs.IdConvert:main",
            "interaction-score = pairs.InteractionScore:main",
            "multi-interaction-score = pairs.MultiInteractionScore:main",
            "random-sequences = pairs.RandomSequences:main",
            "score-matrix = pairs.ScoreMatrix:main",
            "split-fasta = pairs.SplitFasta:main",
        ]
    }
)
