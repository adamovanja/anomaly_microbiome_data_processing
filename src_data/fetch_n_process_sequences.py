import argparse
import subprocess

# from srcd.cluster_sequences import cluster_sequences
# from srcd.denoise_sequences import denoise_sequences
from srcd.trim_sequences import trim_sequences


def _fetch_sequences(version, path_to_data):
    command = (
        f"srcd/fetch_sequences.sh {path_to_data}/runids_v{version} "
        f"{path_to_data}  6"
    )
    subprocess.run(command, shell=True)


def fetch_n_process_sequences(version, path_to_data: str = "../data/raw_data"):
    # fetch
    _fetch_sequences(version, path_to_data)

    # trim
    path2md = f"{path_to_data}/metadata_proc_v{version}.tsv"
    trim_sequences(path2md=path2md, path2seq=path_to_data, threads=6)

    # # denoise
    # path_trunc_len = f"{path_to_data}/trunc_len.json"
    # denoise_sequences(
    #     path2md=path2md, path2trunc_len=path_trunc_len, path2seq=path_to_data, threads=6
    # )

    # # cluster
    # cluster_sequences(path2md=path2md, path2seq=path_to_data, threads=6)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch and process sequences.")
    parser.add_argument(
        "--tag", type=str, required=True, help="Tag for fetching sequences."
    )
    args = parser.parse_args()
    fetch_n_process_sequences(args.tag)
