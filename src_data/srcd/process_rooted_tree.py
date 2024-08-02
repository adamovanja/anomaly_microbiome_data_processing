import argparse
import os

import numpy as np
import skbio
from q2_types.feature_data import DNAIterator
from qiime2 import Artifact
import qiime2 as q2
import pandas as pd 
import skbio
import re

def process_rooted_tree(tree_file):
    phylogeny = q2.Artifact.load(tree_file)

    tree = phylogeny.view(skbio.TreeNode)

    # rename - removing tax info from leaves
    for node in tree.tips():
        # Extract the desired part of the node name using a regular expression
        match = re.match(r"([A-Za-z0-9_\.]+)\s+.*", node.name)
        if match:
            # Update the node name with the extracted part
            node.name = match.group(1)
    
    phylogeny_renamed = Artifact.import_data("Phylogeny[Rooted]", tree)
    file_out = tree_file.replace("rooted-tree", "rooted-tree-proc")
    phylogeny_renamed.save(file_out)
    

if __name__ == "__main__":
    # get user inputs
    parser = argparse.ArgumentParser(description="Remove taxonomic info from nodes.")
    parser.add_argument(
        "--tree",
        type=str,
        required=True,
        help="Location of rooted tree",
    )
    args = parser.parse_args()
    process_rooted_tree(tree_file=args.tree)
