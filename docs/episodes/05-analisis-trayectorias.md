---
layout: default
title: Episode 5 - Trajectory analysis with MDTraj
permalink: /episodes/05-analisis-trayectorias/
---

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/04-muestreo-avanzado/">Previous</a>
  <a href="{{ site.baseurl }}/episodes/">All episodes</a>
  <a href="{{ site.baseurl }}/episodes/06-pyemma/">Next</a>
</div>

<!-- toc:start -->
## Table of contents
- [Duration](#duration)
- [Objectives](#objectives)
- [Content](#content)
- [Alanine dipeptide](#alanine-dipeptide)
- [Protein-ligand complex](#protein-ligand-complex)
- [Simulated annealing](#simulated-annealing)
- [VMD command-line analysis](#vmd-command-line-analysis)
<!-- toc:end -->

## Duration

- **Session:** 60 min
- **Exercises:** 45 min

## Objectives

- Analyze simple and complex trajectories.
- Reimage and align trajectories.
- Generate RMSD with MDTraj.

## Content

- Read DCD files with MDTraj.
- Reimage and align to the initial frame.
- Plots to compare systems.

## Alanine dipeptide

### Guided demo

<!-- sync-from: docs/episodes/scripts/04-analisis-trayectorias_simple.py -->
```python
#!/usr/bin/env python3
import argparse
import os
from pathlib import Path

import mdtraj as md
import matplotlib.pyplot as plt
COURSE_DIR = Path(os.environ.get("COURSE_DIR", str(Path.home() / "Concepcion26"))).expanduser()


def main() -> None:
    parser = argparse.ArgumentParser(description="Analyse alanine dipeptide trajectory")
    parser.add_argument("-p", "--protein", default="alanine-dipeptide.pdb", help="PDB file")
    parser.add_argument("-t", "--trajectory", default="traj.dcd", help="Trajectory DCD file")
    parser.add_argument("-o", "--output", default="rmsd_alanine.png", help="Output figure")
    args = parser.parse_args()

    out_dir = COURSE_DIR / "results" / "04-analisis-trayectorias" / "simple"
    out_dir.mkdir(parents=True, exist_ok=True)
    output_path = out_dir / args.output

    traj = md.load(args.trajectory, top=args.protein)
    rmsd = md.rmsd(traj, traj, frame=0)

    plt.plot(rmsd)
    plt.xlabel("Frame")
    plt.ylabel("RMSD (nm)")
    plt.title("Alanine dipeptide RMSD")
    plt.savefig(output_path, dpi=150)
    print("Written", output_path)


if __name__ == "__main__":
    main()
```

Script source: <a href="{{ site.baseurl }}/episodes/scripts/04-analisis-trayectorias_simple.py">04-analisis-trayectorias_simple.py</a>

Script source: <a href="{{ site.baseurl }}/episodes/scripts/04-analisis-trayectorias_simple.py">04-analisis-trayectorias_simple.py</a>

Script source: <a href="{{ site.baseurl }}/episodes/scripts/04-analisis-trayectorias_simple.py">04-analisis-trayectorias_simple.py</a>

Script source: <a href="{{ site.baseurl }}/episodes/scripts/04-analisis-trayectorias_simple.py">04-analisis-trayectorias_simple.py</a>

Script source: <a href="{{ site.baseurl }}/episodes/scripts/04-analisis-trayectorias_simple.py">04-analisis-trayectorias_simple.py</a>

Script source: <a href="{{ site.baseurl }}/episodes/scripts/04-analisis-trayectorias_simple.py">04-analisis-trayectorias_simple.py</a>

### Exercise

- Run the simple analysis with <a href="{{ site.baseurl }}/data/">traj.dcd</a>.

### Key points

- The metrics are comparable across systems.

### Notebooks and scripts

- This notebook reads the alanine trajectory with MDTraj, reimages frames, and compares RMSD metrics for the simple system. (<a href="{{ site.baseurl }}/episodes/notebooks/05-analisis-trayectorias_simple.ipynb">notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/04-analisis-trayectorias_simple.py">script</a>)

## Protein-ligand complex

### Guided demo

<!-- sync-from: docs/episodes/scripts/04-analisis-trayectorias.py -->
```python
#!/usr/bin/env python3
import argparse
import os
from pathlib import Path

import mdtraj as md
import plotly.graph_objects as go
COURSE_DIR = Path(os.environ.get("COURSE_DIR", str(Path.home() / "Concepcion26"))).expanduser()


def main() -> None:
    parser = argparse.ArgumentParser(description="Analyse complex trajectory")
    parser.add_argument("-p", "--protein", default="output_minimised.pdb", help="Protein PDB file")
    parser.add_argument("-t", "--trajectory", default="output_traj.dcd", help="Trajectory DCD file")
    parser.add_argument("-o", "--output", default="output_reimaged", help="Output base name")
    parser.add_argument("-r", "--remove-waters", action="store_true", help="Remove waters and ions")
    args = parser.parse_args(args=[])

    print("Reading trajectory", args.trajectory)
    traj = md.load(args.trajectory, top=args.protein)
    traj.image_molecules(inplace=True)

    if args.remove_waters:
        print("Removing waters")
        traj = traj.atom_slice(traj.top.select("not resname HOH POPC CL NA"))

    print("Realigning")
    prot = traj.top.select("protein")
    traj.superpose(traj[0], atom_indices=prot)

    out_dir = COURSE_DIR / "results" / "04-analisis-trayectorias" / "complex"
    out_dir.mkdir(parents=True, exist_ok=True)
    output_base = str(out_dir / args.output)

    print("Writing re-imaged PDB", output_base + ".pdb")
    traj[0].save(output_base + ".pdb")

    print("Writing re-imaged trajectory", output_base + ".dcd")
    traj.save(output_base + ".dcd")

    atoms = traj.top.select("chainid 1")
    rmsd_lig = md.rmsd(traj, traj, frame=0, atom_indices=atoms, parallel=True, precentered=False)

    atoms = traj.top.select("chainid 0 and backbone")
    rmsd_bck = md.rmsd(traj, traj, frame=0, atom_indices=atoms, parallel=True, precentered=False)

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=traj.time, y=rmsd_lig, mode="lines", name="Ligand"))
    fig.add_trace(go.Scatter(x=traj.time, y=rmsd_bck, mode="lines", name="Backbone"))

    fig.update_layout(title="Trajectory for " + args.trajectory, xaxis_title="Frame", yaxis_title="RMSD")

    file = output_base + ".svg"
    print("Writing RMSD output to", file)
    fig.write_image(file)


if __name__ == "__main__":
    main()
```

Script source: <a href="{{ site.baseurl }}/episodes/scripts/04-analisis-trayectorias.py">04-analisis-trayectorias.py</a>

Script source: <a href="{{ site.baseurl }}/episodes/scripts/04-analisis-trayectorias.py">04-analisis-trayectorias.py</a>

Script source: <a href="{{ site.baseurl }}/episodes/scripts/04-analisis-trayectorias.py">04-analisis-trayectorias.py</a>

Script source: <a href="{{ site.baseurl }}/episodes/scripts/04-analisis-trayectorias.py">04-analisis-trayectorias.py</a>

Script source: <a href="{{ site.baseurl }}/episodes/scripts/04-analisis-trayectorias.py">04-analisis-trayectorias.py</a>

Script source: <a href="{{ site.baseurl }}/episodes/scripts/04-analisis-trayectorias.py">04-analisis-trayectorias.py</a>

### Exercise

- Run the complex analysis with <a href="{{ site.baseurl }}/reference/">output_traj.dcd</a>.
- Save the resulting figures.

### Key points

- MDTraj enables reproducible analysis with a few lines.

### Notebooks and scripts

- This notebook performs MDTraj analyses on the protein-ligand trajectory, re-imaging the complex and exporting RMSD plots to highlight differences. (<a href="{{ site.baseurl }}/episodes/notebooks/05-analisis-trayectorias.ipynb">notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/04-analisis-trayectorias.py">script</a>)

## Simulated annealing

### Guided demo

<!-- sync-from: docs/episodes/notebooks/04-simulated-annealing.ipynb -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/04-simulated-annealing.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/04-simulated-annealing.ipynb" download>Download notebook</a></div></div>

### Exercise

- Run the notebook from `COURSE_DIR` so the annealing loop reads `data/alanine-dipeptide.pdb` and writes snapshots under `results/04-annealing/`; tweak `steps_per_temp` or the number of increments to sample more aggressive schedules.

### Key points

- Simulated annealing also lives inside the `COURSE_DIR` tree, so the notebook and its outputs are consistent with the other Episode 04 materials.

### Notebooks

- Simulated annealing with OpenMM plus the course-standard datasets. (<a href="{{ site.baseurl }}/episodes/notebooks/04-simulated-annealing.ipynb">notebook</a>)

## VMD command-line analysis

We also provide a lightweight VMD Tcl helper that runs from the terminal to compute RMSD, RMSF, specific distances, and even export the movie frames for the complex trajectory. The script lives at `docs/episodes/scripts/04-vmd-analysis.tcl` and requires `$COURSE_DIR/results/03-simulaciones-clasicas/complex/output_traj.dcd`.

Execute it with:

```
COURSE_DIR=~/Concepcion26 vmd -dispdev text -e docs/episodes/scripts/04-vmd-analysis.tcl
```

The run writes CSV summaries into `$COURSE_DIR/results/04-analisis-trayectorias/complex/vmd` and emits a PNG sequence under `movie/` (Tachyon renders each frame). The distance file tracks the anchor residue/ligand separation frame by frame, giving a quick consistency check.

If you’d like a guided introduction before using the script, review the VMD “Analysis of MD Trajectories” tutorial at https://www.ks.uiuc.edu/Training/Tutorials/vmd/analysis/ — the commands there are the same ones wrapped inside our Tcl helper, and the pages explain how to instrument RMSD/RMSF/distance traces interactively.

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/04-muestreo-avanzado/">Previous</a>
  <a href="{{ site.baseurl }}/episodes/">All episodes</a>
  <a href="{{ site.baseurl }}/episodes/06-pyemma/">Next</a>
</div>
