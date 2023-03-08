# RCA

**RCA Pipeline**

RCA is a pipeline for the discovery of viruses from Illumina and Nanopore rolling circle amplification data. 

**Setting up a Conda Environment**

The RCA.yml file provided allows for easily setting up a conda environment, which will contain most of the tools needed to run the script.
Ensure that conda is installed and active, and then execute `conda env create -f RCA.yml` to create and install all tools into the environment.

Please install Ashure separately and update the `--ashure_path` parameter with the location of the ashure.py file.

Please download the blast database in the ICTV folder and point the `--ictv_db` parameter to the ICTV_viruses.fasta file.

Please activate the conda environment before running the scripts.

Usage:

`python Scripts/RCA/rca.py illumina --input Illumina_run/raw_data/ --output Illumina_run/`

`python Scripts/RCA/rca.py nanopore --input Nanopore_run/raw_data/ --output Nanopore_run/`

**Illumina Parameters**
```
usage: RCA illumina [-h] --input INPUT --output OUTPUT [--create_dirs {Y,N}]
                    [--bbduk {Y,N}] [--spades {Y,N}] [--virus_detection {Y,N}]
                    [--doc_env {Y,N}] [--adapters ADAPTERS]
                    [--min_len MIN_LEN] [--min_quality MIN_QUALITY]
                    [--threads THREADS] [--memory MEMORY] [--hmm HMM]

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT         This is the location of the raw data directory.
  --output OUTPUT       This is where the output data will be generated.
  --create_dirs {Y,N}   Creates the directory structure. Default Y.
  --bbduk {Y,N}         Runs Bbduk. Default Y.
  --spades {Y,N}        Runs isONclust. Default Y.
  --virus_detection {Y,N}
                        Extract potential viruses. Default Y.
  --doc_env {Y,N}       Creates environment and parameter details. Default Y.
  --adapters ADAPTERS   bbduk adapter references.
  --min_len MIN_LEN     This is the mininum length of a read allowed by bbduk.
                        Default 50.
  --min_quality MIN_QUALITY
                        bbduk phred quality trim parameter. Default 10.
  --threads THREADS     The number of threads to run Spades with. Default 60.
  --memory MEMORY       The amount of memory, in GB, to run Spades with.
                        Default 500.
  --hmm HMM             This is the path to the hmm file.
  ```
  
  **Nanopore Parameters**
  ```
  usage: RCA nanopore [-h] --input INPUT --output OUTPUT [--ictv_db ICTV_DB]
                    [--ashure_path ASHURE_PATH] [--create_dirs {Y,N}]
                    [--ashure {Y,N}] [--isonclust {Y,N}] [--medaka {Y,N}]
                    [--doc_env {Y,N}] [--threads THREADS]
                    [--quality_score QUALITY_SCORE]

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT         This is the location of the raw data directory.
  --output OUTPUT       This is where the output data will be generated.
  --ictv_db ICTV_DB     This is the path to the ictv database fasta file.
  --ashure_path ASHURE_PATH
                        This is the path to the ashure script.
  --create_dirs {Y,N}   Creates the directory structure. Default Y.
  --ashure {Y,N}        Runs Ashure. Default Y.
  --isonclust {Y,N}     Runs isONclust. Default Y.
  --medaka {Y,N}        Runs medaka. Default Y.
  --doc_env {Y,N}       Creates environment and parameter details. Default Y.
  --threads THREADS     The number of threads to run isONclust and Medaka
                        with. Default 120.
  --quality_score QUALITY_SCORE
                        The quality score cutoff for isONclust. Default 7.
```
