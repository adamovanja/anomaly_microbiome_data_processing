import argparse
import subprocess

from src.cluster_sequences import cluster_sequences
from src.denoise_sequences import denoise_sequences
from src.trim_sequences import trim_sequences
from src.check_sequences import check_sequences


def _fetch_sequences(n_threads, path_to_data):
    command = (
        f"../src/fetch_sequences.sh {path_to_data}/runids"
        f"{path_to_data}  {n_threads}"
    )
    subprocess.run(command, shell=True)


def fetch_n_process_sequences(version, n_threads, path_to_data: str = "../data/raw"):
    # fetch
    _fetch_sequences(n_threads, path_to_data)

    # trim
    path2md = f"{path_to_data}/metadata_proc_v{version}.tsv"
    trim_sequences(path2md=path2md, path2seq=path_to_data, threads=n_threads)

    # denoise
    path_trunc_len = f"{path_to_data}/trunc_len.json"
    denoise_sequences(
        path2md=path2md,
        path2trunc_len=path_trunc_len,
        path2seq=path_to_data,
        threads=n_threads,
    )

    # cluster
    cluster_sequences(path2md=path2md, path2seq=path_to_data, threads=n_threads)

    # save stats of all processed sequences - can be viewed in
    # describe_sequences.ipynb
    check_sequences(path2md=path2md)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch and process sequences.")
    parser.add_argument(
        "--tag", type=str, required=True, help="Tag for fetching sequences."
    )
    parser.add_argument(
        "--n_threads", type=int, required=True, help="Number of threads to use."
    )
    args = parser.parse_args()
    fetch_n_process_sequences(args.tag, args.n_threads)
