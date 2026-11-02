"""Tests for the bootstrapped alpha diversity of rarefied feature tables."""

import numpy as np
import pandas as pd
import qiime2 as q2
import skbio

from src import postprocess_sequences as ps

TREE_FILE = "silva-138-1-ssu-nr99-seqs-rooted-tree-proc-ben.qza"


def _feature_table(n_samples=3, n_features=40, reads=4000, seed=0):
    rng = np.random.default_rng(seed)
    counts = rng.multinomial(reads, rng.dirichlet(np.ones(n_features)), size=n_samples)
    return pd.DataFrame(
        counts.astype(float),
        index=[f"s{i}" for i in range(n_samples)],
        columns=[f"f{j}" for j in range(n_features)],
    )


def test_worker_count_leaves_two_cpus_free(monkeypatch):
    monkeypatch.setattr(ps.os, "cpu_count", lambda: 8)
    assert ps._n_workers(10, 500) == 6
    assert ps._n_workers(4, 500) == 4
    assert ps._n_workers(10, 3) == 3
    monkeypatch.setattr(ps.os, "cpu_count", lambda: 2)
    assert ps._n_workers(10, 500) == 1


def test_rarefied_alpha_runs_are_independent_rarefactions_across_workers():
    runs = ps.rarefied_alpha_runs(
        "shannon", _feature_table(), 100, None, n_threads=2, n=6
    )

    assert runs.index.tolist() == ["s0", "s1", "s2"]
    assert runs.columns.tolist() == [f"shannon_{i}" for i in range(6)]
    # no two rarefactions are identical, also not those of different workers
    assert runs.T.drop_duplicates().shape[0] == 6


def test_bootstrapped_alpha_div_averages_the_rarefactions():
    alpha = ps.bootstrapped_alpha_div(
        "observed_features", _feature_table(), 100, None, 2, n=4
    )

    assert alpha.name == "div_alpha_observed_features"
    assert alpha.index.tolist() == ["s0", "s1", "s2"]
    assert ((alpha > 0) & (alpha <= 40)).all()


def test_faith_pd_runs_use_the_tree_in_path(tmp_path):
    # four equally abundant features: all are kept when rarefying 100 of
    # 4000 reads, so every run spans the whole tree (4 tips + 2 inner edges)
    ft = pd.DataFrame(1000.0, index=["s0", "s1"], columns=["f0", "f1", "f2", "f3"])
    tree = skbio.TreeNode.read(["((f0:1,f1:1):1,(f2:1,f3:1):1);"])
    q2.Artifact.import_data("Phylogeny[Rooted]", tree).save(str(tmp_path / TREE_FILE))

    runs = ps.rarefied_alpha_runs("faith_pd", ft, 100, str(tmp_path), n_threads=2, n=4)

    assert runs.shape == (2, 4)
    assert (runs == 6.0).all().all()
