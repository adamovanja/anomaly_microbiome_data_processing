# anomaly_microbiome_data_processing
This repos processes the microbial feature tables as needed for the anomaly detection usecase.

## Setup
To fetch and preprocess the metadata and sequences of the cohort, we suggest you set up a conda environment as follows:
* Install [mamba](https://github.com/mamba-org/mamba)
* Create and activate a conda environment with the required dependencies:
```shell
mamba env create -f environment.yml
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
cd scripts
python fetch_n_process_metadata.py --email your@mail.com --n_threads 20
````
Beware: When running this command for the first time all metadata and the raw amplicon nucleotide sequences need to be fetched from NCBI SRA (~5-20 min). Due to the restrictions of some journal's website, not all of the required supplementary material can be fetched programmatically. In this case, follow the ValueError instructions and re-execute the above command.

### Get sequences
To fetch and process the respective sequences insert the tag outputted at the end of the previous metadata script as `$TAG`, select a number of threads `$THREADS` and run:

````
cd scripts
python fetch_n_process_sequences.py --tag $TAG --n_threads $THREADS
````
Beware: Fetching and processing raw nucleotide sequences requires a total of ~70 GB storage space. The processing time depends on your internet connection speed and your computing power. As a reference, here is the approximate duration for each step when performed with `n_threads 6` on a MacOS system equipped with a 2 GHz Quad-Core Intel Core i5 processor and 16 GB of 3733 MHz LPDDR4X memory:
fetching: 24 hrs
trimming: 2 hrs
denoising: 6 hrs
clustering: 13 hrs

And on an Ubuntu-based high-performance cluster with 200GB RAM and 50 CPU cores selected:
fetching: 24 hrs
trimming: x hrs
denoising: x hrs
clustering: x hrs

### Create feature table used for modelling
Once the previous commands worked successfully you can create the final feature table to be used for modelling using the same `$TAG` you used before. When running this command for the first time a phylogenetic tree is build from the aligned SILVA reference sequence database, which requires ~200GB of memory:
````
cd scripts
python create_feature_table.py --tag $TAG  --n_threads $THREADS
````

### Down-stream data analyses

#### Description of the resulting feature table
... can be created in the notebook `scripts/describe_dataset.ipynb`.

#### Analysis of matched alpha diversity after antibiotics exposure
... can be created in the notebook `scripts/describe_matched_alpha.ipynb`.

#### Description of sequencing length and sample counts
... after each sequence processing step can be created in the notebook `scripts/describe_sequences.ipynb`.
