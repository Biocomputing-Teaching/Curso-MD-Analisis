from pathlib import Path
import os

from contextlib import nullcontext
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import mdshare
import pyemma

COURSE_DIR = Path(os.environ.get('COURSE_DIR', str(Path.home() / 'Concepcion26'))).expanduser()
DATA_DIR = COURSE_DIR / 'data'
DATA_DIR.mkdir(parents=True, exist_ok=True)

# %%

pdb = mdshare.fetch('pentapeptide-impl-solv.pdb', working_directory=str(DATA_DIR))
files = mdshare.fetch('pentapeptide-*-500ns-impl-solv.xtc', working_directory=str(DATA_DIR))
print('Topology:', pdb)
print('Trajectories:', len(files))

# %%

torsions_feat = pyemma.coordinates.featurizer(str(pdb))
torsions_feat.add_backbone_torsions(cossin=True, periodic=False)
torsions_data = pyemma.coordinates.load(files, features=torsions_feat)

positions_feat = pyemma.coordinates.featurizer(str(pdb))
positions_feat.add_selection(positions_feat.select_Backbone())
positions_data = pyemma.coordinates.load(files, features=positions_feat)

backbone = torsions_feat.select_Backbone()
distances_feat = pyemma.coordinates.featurizer(str(pdb))
distances_feat.add_distances(distances_feat.pairs(backbone, excluded_neighbors=2), periodic=False)
distances_data = pyemma.coordinates.load(files, features=distances_feat)

# %%

def score_cv(data, dim, lag, number_of_splits=10, validation_fraction=0.5):
    contexts = getattr(pyemma.util, 'contexts', None)
    if contexts is None:
        ctx = nullcontext()
    else:
        ctx = contexts.settings(show_progress_bars=False)
    with ctx:
        nval = int(len(data) * validation_fraction)
        scores = np.zeros(number_of_splits)
        for n in range(number_of_splits):
            ival = np.random.choice(len(data), size=nval, replace=False)
            vamp = pyemma.coordinates.vamp([
                d for i, d in enumerate(data) if i not in ival], lag=lag, dim=dim)
            scores[n] = vamp.score([d for i, d in enumerate(data) if i in ival])
    return scores

labels = ['torsions', 'positions', 'distances']

fig, axes = plt.subplots(1, 3, figsize=(12, 3), sharey=True)
for ax, lag in zip(axes.flat, [5, 10, 20]):
    data_scores = []
    data_errors = []
    for dataset in (torsions_data, positions_data, distances_data):
        sc = score_cv(dataset, dim=10, lag=lag)
        data_scores.append(sc.mean())
        data_errors.append(sc.std())
    ax.bar(labels, data_scores, yerr=data_errors)
    ax.set_title(f'lag {lag * 0.1:.1f} ns')
axes[0].set_ylabel('VAMP2 score')
fig.tight_layout()

# %%

tica = pyemma.coordinates.tica(torsions_data, lag=5)
tica_output = tica.get_output()
tica_concatenated = np.concatenate(tica_output)

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
pyemma.plots.plot_feature_histograms(tica_concatenated, ax=axes[0], feature_labels=[f'IC{i+1}' for i in range(4)], ylog=True)
pyemma.plots.plot_density(*tica_concatenated[:, :2].T, ax=axes[1], logscale=True)
axes[1].set_xlabel('IC 1')
axes[1].set_ylabel('IC 2')
fig.tight_layout()

# %%

ncenters = [5, 10, 30, 75, 200, 450]
scores = np.zeros((len(ncenters), 5))
contexts = getattr(pyemma.util, 'contexts', None)
for idx, k in enumerate(ncenters):
    for trial in range(5):
        ctx = contexts.settings(show_progress_bars=False) if contexts else nullcontext()
        with ctx:
            cl = pyemma.coordinates.cluster_kmeans(tica_output, k=k, max_iter=50, stride=50)
            msm = pyemma.msm.estimate_markov_model(cl.dtrajs, 5)
            scores[idx, trial] = msm.score_cv(cl.dtrajs, n=1, score_method='VAMP2', score_k=min(10, k))

fig, ax = plt.subplots()
lower, upper = pyemma.util.statistics.confidence_interval(scores.T.tolist(), conf=0.9)
ax.fill_between(ncenters, lower, upper, alpha=0.3)
ax.plot(ncenters, scores.mean(axis=1), '-o')
ax.set_xscale('log')
ax.set_xlabel('number of cluster centers')
ax.set_ylabel('VAMP-2 score')
fig.tight_layout()

k = 75
cluster = pyemma.coordinates.cluster_kmeans(tica_output, k=k, max_iter=50, stride=10, fixed_seed=1)
dtrajs = np.concatenate(cluster.dtrajs)

# %%

msm = pyemma.msm.bayesian_markov_model(cluster.dtrajs, lag=5, dt_traj='0.1 ns')
print('fraction of states used:', msm.active_state_fraction)
print('fraction of counts used:', msm.active_count_fraction)
its = pyemma.msm.its(cluster.dtrajs, lags=50, nits=15, errors='bayes')
pyemma.plots.plot_implied_timescales(its, units='ns', dt=0.1)

# %%

nstates = 5
cktest = msm.cktest(nstates, mlags=6)
pyemma.plots.plot_cktest(cktest, dt=0.1, units='ns')

# %%

eigvec = msm.eigenvectors_right()
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
pyemma.plots.plot_contour(*tica_concatenated[:, :2].T, eigvec[dtrajs, 1], ax=axes[0], cmap='PiYG', mask=True)
pyemma.plots.plot_contour(*tica_concatenated[:, :2].T, eigvec[dtrajs, 2], ax=axes[1], cmap='PiYG', mask=True)
[ax.set_xlabel('IC 1') for ax in axes];
axes[0].set_ylabel('IC 2')

# %%

pcca = msm.pcca(nstates)
metatraj = msm.metastable_assignments[dtrajs]
fig, ax = plt.subplots(figsize=(4, 4))
pyemma.plots.plot_state_map(*tica_concatenated[:, :2].T, metatraj, ax=ax)
ax.set_xlabel('IC 1')
ax.set_ylabel('IC 2')

# %%

start, final = 1, 3
A = msm.metastable_sets[start]
B = msm.metastable_sets[final]
flux = pyemma.msm.tpt(msm, A, B)
cg, cgflux = flux.coarse_grain(msm.metastable_sets)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3.5))
pyemma.plots.plot_contour(
    *tica_concatenated[:, :2].T,
    flux.committor[dtrajs],
    ax=ax1,
    mask=True,
    cmap="brg",
    cbar_label="committor $\mathcal{S}_2 \to \mathcal{S}_4$")
ax1.set_xlabel("IC 1")
ax1.set_ylabel("IC 2")
pyemma.plots.plot_flux(
    cgflux,
    ax=ax2,
    show_committor=True,
    max_width=15,
    max_height=15,
    minflux=1e-5)
ax2.set_yticks([])
ax2.set_xlabel("coarse-grained states")
fig.tight_layout()

# %%

from itertools import product
from mdtraj import shrake_rupley, compute_rg

markov_samples = [smpl for smpl in msm.sample_by_state(20)]
reader = pyemma.coordinates.source(files, top=str(pdb))
samples = [pyemma.coordinates.save_traj(reader, sample, outfile=None, top=str(pdb)) for sample in markov_samples]
markov_sasa = [shrake_rupley(sample, mode='residue') for sample in samples]
markov_rg = [compute_rg(sample) for sample in samples]
markov_average_rg = np.array(markov_rg).mean(axis=1)
print('Average Rg (nm):', msm.expectation(markov_average_rg))
