#!/usr/bin/env python3
import argparse
import os
from pathlib import Path
import numpy as np
import pyemma
import matplotlib.pyplot as plt
import pyemma.plots


def load_features(topology: str):
    feat = pyemma.coordinates.featurizer(topology)
    feat.add_backbone_torsions(periodic=False)
    return feat


def run_pipeline(traj_files, topology, lag, k):
    feat = load_features(topology)
    data = pyemma.coordinates.load(traj_files, features=feat)

    tica = pyemma.coordinates.tica(data, lag=lag)
    tica_output = tica.get_output()

    cluster = pyemma.coordinates.cluster_kmeans(tica_output, k=k, max_iter=50)
    msm = pyemma.msm.estimate_markov_model(cluster.dtrajs, lag=lag)

    return data, tica, cluster, msm


def main() -> None:
    parser = argparse.ArgumentParser(description="PyEMMA tutorials on alanine dipeptide")
    parser.add_argument("-t", "--traj", nargs="+", default=None, help="Trajectory files")
    parser.add_argument("-p", "--top", default=None, help="Topology PDB")
    parser.add_argument("--lag", type=int, default=10, help="Lag time")
    parser.add_argument("-k", type=int, default=50, help="Number of clusters")
    args = parser.parse_args()

    data_dir = COURSE_DIR / "data"
    traj_files = args.traj or [str(data_dir / "alanine-dipeptide.dcd") ]
    topology = args.top or str(data_dir / "alanine-dipeptide.pdb")

    out_dir = COURSE_DIR / "results" / "06-pyemma" / "simple"
    out_dir.mkdir(parents=True, exist_ok=True)
    data, tica, cluster, msm = run_pipeline(traj_files, topology, args.lag, args.k)

    print("Tutorial 01 - Data IO & featurization: OK")
    print("Feature dimension:", tica.dimension())

    print("Tutorial 02 - TICA & clustering: OK")
    print("TICA lag:", tica.lag)

    print("Tutorial 03 - MSM estimation & validation")
    its = pyemma.msm.its(cluster.dtrajs, lags=[1, 2, 5, 10])
    print("ITS computed for", len(its.lags), "lags")
    plt.figure()
    pyemma.plots.plot_implied_timescales(its)
    plt.savefig(out_dir / "its.png", dpi=150)

    print("Tutorial 04 - MSM analysis")
    timescales = msm.timescales()[:5]
    print("Timescales:", timescales)
    (out_dir / "timescales.txt").write_text("\\n".join(str(t) for t in timescales))

    print("Tutorial 05 - PCCA & TPT")
    try:
        msm.pcca(2)
        tpt = pyemma.msm.tpt(msm, msm.active_set[0], msm.active_set[-1])
        print("TPT flux:", tpt.total_flux)
    except Exception as exc:
        print("PCCA/TPT skipped:", exc)

    print("Tutorial 06 - Expectations & observables")
    try:
        obs = data[0][:, 0]
        print("Observable mean:", float(np.mean(obs)))
    except Exception as exc:
        print("Observables skipped:", exc)

    print("Tutorial 07 - HMM")
    try:
        hmm = pyemma.msm.bayesian_hidden_markov_model(cluster.dtrajs, nstates=2, lag=args.lag)
        print("HMM timescales:", hmm.timescales()[:3])
    except Exception as exc:
        print("HMM skipped:", exc)

    print("Tutorial 08 - Common problems")
    print("Active fraction:", msm.active_state_fraction)
    print("Connected count fraction:", msm.active_count_fraction)


if __name__ == "__main__":
    main()
