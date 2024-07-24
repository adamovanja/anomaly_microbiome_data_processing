# anomaly_microbiome_data_processing


## Setup
To fetch and preprocess the metadata and sequences of the cohort, we suggest you set up a conda environment as follows:
* Install [mamba](https://github.com/mamba-org/mamba)
* Create and activate a conda environment with the required dependencies:
```shell
cd src_data
mamba env create -f data_environment.yml
conda activate preprocess_microbiome
conda install -c https://packages.qiime2.org/qiime2/2024.5/metagenome/released/ -c conda-forge -c bioconda -c defaults q2-fondue -y
pip install -e .
qiime dev refresh-cache
```
* Run the `vdb-config` tool and exit by pressing x (needed to initialize the wrapped SRA Toolkit for more information see [here](https://github.com/ncbi/sra-tools/wiki/05.-Toolkit-Configuration))
```shell
vdb-config -i
```

* In case you need to configure a proxy server, run the following command:
```shell
vdb-config --proxy <your proxy URL> --proxy-disable no
```
## Creating the processed sequences and metadata 
Run in your terminal:
````
python fetch_n_process_metadata.py
````
Beware: When running this command for the first time all metadata and the raw amplicon nucleotide sequences need to be fetched from NCBI SRA. This takes roughly 20min on a MacOS with `n_jobs=6`. Also, due to the restrictions of some journal's website, not all of the required supplementary material can be fetched programmatically. In this case, follow the ValueError instructions and re-execute the above command.
