#!/usr/bin/env python3
import argparse
import os
from pathlib import Path

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from deeptime.clustering import KMeans
from deeptime.markov.msm import MaximumLikelihoodMSM

def compute_features(topology: str, trajectories: list[str]) -> np.ndarray:
    traj = md.load(trajectories, top=topology)
    atom_indices = traj.top.select("backbone")
    pairs = [(atom_indices[i], atom_indices[i + 1]) for i in range(len(atom_indices) - 1)]
    distances = md.compute_distances(traj, pairs)
    return distances


def main() -> None:
    parser = argparse.ArgumentParser(description="deeptime MSM on protein-ligand complex")
    parser.add_argument("-t", "--traj", nargs="+", default=None, help="Trajectory files")
    parser.add_argument("-p", "--top", default=None, help="Topology file")
    parser.add_argument("-k", "--k", type=int, default=50, help="Number of clusters")
    parser.add_argument("--lag", type=int, default=10, help="Lag time")
    parser.add_argument("--pca", type=int, default=2, help="PCA components")
    args = parser.parse_args()

    data_dir = COURSE_DIR / "data" / "complex"
    traj_files = args.traj or [str(data_dir / "output_traj.dcd") ]
    topology = args.top or str(data_dir / "output_minimised.pdb")

    out_dir = COURSE_DIR / "results" / "07-deeptime" / "complex"
    out_dir.mkdir(parents=True, exist_ok=True)
    features = compute_features(topology, traj_files)
    pca = PCA(n_components=args.pca)
    pca_coords = pca.fit_transform(features)

    plt.figure()
    plt.scatter(pca_coords[:, 0], pca_coords[:, 1], s=5, alpha=0.5)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title("Complex PCA")
    plt.savefig(out_dir / "pca.png", dpi=150)

    kmeans = KMeans(n_clusters=args.k, max_iter=50)
    dtrajs = kmeans.fit_transform(pca_coords)

    msm = MaximumLikelihoodMSM(lagtime=args.lag).fit(dtrajs).fetch_model()
    timescales = msm.timescales()[:5]
    print("Timescales:", timescales)
    (out_dir / "timescales.txt").write_text("\\n".join(str(t) for t in timescales))


if __name__ == "__main__":
    main()
