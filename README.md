# anomaly_microbiome_data_processing
This repos processes the microbial feature tables as needed for the anomaly detection usecase.

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

### Get metadata
To fetch and process the metadata run:
````
cd src_data
python fetch_n_process_metadata.py --email your@mail.com --n_threads 20
````
Beware: When running this command for the first time all metadata and the raw amplicon nucleotide sequences need to be fetched from NCBI SRA (~5-20 min). Due to the restrictions of some journal's website, not all of the required supplementary material can be fetched programmatically. In this case, follow the ValueError instructions and re-execute the above command.

### Get sequences
To fetch and process the respective sequences insert the tag outputted at the end of the previous metadata script as `$TAG` and run:

````
cd src_data
python fetch_n_process_sequences.py --tag $TAG --n_threads 50
````
Beware: Fetching and processing raw nucleotide sequences requires ~40 GB of storage space. The processing time depends on your internet connection speed and your computing power - as a reference here is the approximate duration for each step performed on HPC with ~27s/it, 50 CPUs (Memory per CPU=4096MB) and 50 threads selected:

fetching: 27 hrs
trimming: 1 hr
denoising: tba
clustering: tba

### Create feature table used for modelling
Once the previous commands worked successfully you can create the final feature table to be used for modelling using the same `$TAG` you used before:
````
cd src_data
python create_feature_table.py --tag $TAG
````

### Description of the resulting feature table
... can be found in the notebook `src_data/describe_dataset.ipynb`.
