# RCA

**RCA Pipeline**

RCA is a pipeline for the discovery of viruses from Illumina and Nanopore rolling circle amplification data. 

**Setting up a Conda Environment**

The RCA.yml file provided allows for easily setting up a conda environment, which will contain most of the tools needed to run the script.
Ensure that conda is installed and active, and then execute `conda env create -f RCA.yml` to create and install all tools into the environment.

Please install Ashure separately and update the `--ashure_path` parameter with the location of the ashure.py file.

Please activate the conda environment before running the scripts.

Usage:

`python Scripts/RCA/rca.py illumina --input Illumina_run/raw_data/ --output Illumina_run/`
`python Scripts/RCA/rca.py nanopore --input Nanopore_run/raw_data/ --output Nanopore_run/`
