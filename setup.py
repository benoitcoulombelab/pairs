from setuptools import setup, find_packages

setup(
    name="DeepProteinConnector",
    version="1.0-SNAPSHOT",
    packages=find_packages(),
    author="Christian Poitras",
    author_email="christian.poitras@ircm.qc.ca",
    description="Find interactions using AlphaFold",
    keywords="bioinformatics, AlphaFold",
    url="https://github.com/benoitcoulombelab/DeepProteinConnector.git",
    license="GNU General Public License version 3",
    classifiers=[
        "License :: OSI Approved :: GNU General Public License version 3"
    ],
    install_requires=[
        "biopython>=1.81"
    ],
    entry_points={
        "console_scripts": [
            "count-interactions = DeepProteinConnector.CountInteractions:main",
            "fasta-id = DeepProteinConnector.FastaId:main",
            "keep-gene-proteins = DeepProteinConnector.KeepGeneProteins:main",
            "merge-fastas = DeepProteinConnector.MergeFastas:main",
            "split-fasta = DeepProteinConnector.SplitFasta:main"
        ]
    }
)
