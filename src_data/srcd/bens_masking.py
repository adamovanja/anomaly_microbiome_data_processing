import argparse
import os

import numpy as np
import skbio
from q2_types.feature_data import DNAIterator
from qiime2 import Artifact


def bens_masking(file):
    """Perform Ben Kaehler's masking of SILVA aligned reference sequences
    approach from:
    https://gist.github.com/BenKaehler/d9291d59bce5cd3d2a90c73b822b3a21 linked
    here:
    https://forum.qiime2.org/t/estimating-time-and-memory-for-masking-large-dataset/21166/9
    only keeping columns with less than 99.56% gaps.

    Note: this approach is equivalent in result but much faster than `qiime
    alignment mask` with `--p-max-gap-frequency 0.0044`.
    """
    # Read the aligned sequences into a NumPy array for fast access by column.
    silva_aln = Artifact.load(file)
    seq_array = []
    for seq in silva_aln.view(DNAIterator):
        seq_array.append(seq.values)
    seq_array = np.array(seq_array)

    # Count the gaps in each column.
    num_gaps = np.empty(seq_array.shape[1])
    for j in range(seq_array.shape[1]):
        num_gaps[j] = (seq_array[:, j] == b".").sum()
        num_gaps[j] += (seq_array[:, j] == b"-").sum()
        if j % 10000 == 0:
            print(f"Done {j} columns")

    # Find the less gappy columns.
    # predefined number that represents the minimum number of non-gap characters
    # a column should have to be kept:
    num_ok = 1907
    ix = num_gaps <= seq_array.shape[0] - num_ok
    print(
        f"Keeping {ix.sum()} columns that are less than {100 - num_ok/seq_array.shape[0]*100:.2f}% gaps"
    )

    masked = skbio.alignment.TabularMSA(s[ix] for s in silva_aln.view(DNAIterator))
    masked = Artifact.import_data("FeatureData[AlignedSequence]", masked)
    file_out = file.replace(".qza", "-ben-masked.qza")
    masked.save(file_out)


if __name__ == "__main__":
    # get user inputs
    parser = argparse.ArgumentParser(description="Perform Ben Kaehler's masking.")
    parser.add_argument(
        "--seqs",
        type=str,
        required=True,
        help="Location of aligned sequences",
    )
    args = parser.parse_args()
    bens_masking(args.seqs)
