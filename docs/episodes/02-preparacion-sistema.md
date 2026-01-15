---
layout: default
title: Episode 3 - Preparing and editing the system in OpenMM
permalink: /episodes/02-preparacion-sistema/
---

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/03-simulaciones-clasicas/">Previous</a>
  <a href="{{ site.baseurl }}/episodes/">All episodes</a>
  <a href="{{ site.baseurl }}/episodes/05-muestreo-avanzado/">Next</a>
</div>

> **System focus** Episode 3 works with the canonical `protein2-ligand2.pdb` complex (`$COURSE_DIR/data/complex`). Re-run `bash docs/episodes/scripts/01-split-protein2-ligand2.sh` whenever the data directory is refreshed to regenerate `protein2.pdb` and `ligand2.pdb` for every script below.

<!-- toc:start -->
## Table of contents
- [Duration](#duration)
- [Objectives](#objectives)
- [Content](#content)
- [Prepare the working folder](#prepare-the-working-folder)
- [Alanine dipeptide](#alanine-dipeptide)
- [Protein-ligand complex](#protein-ligand-complex)
- [Key points](#key-points)
- [Modeller (OpenMM guide)](#modeller-openmm-guide)
<!-- toc:end -->

## Duration

- **Session:** 60 min
- **Exercises:** 45 min

## Objectives

- Prepare a simple system (alanine) and a complex system (protein-ligand).
- Reproduce the Modeller operations from the OpenMM guide.
- Generate simulation-ready structures and save them reproducibly.

## Content

- Simple preparation with Modeller and light solvation.
- PDBFixer to repair complex proteins.
- Modeller: hydrogens, solvent, membranes, extra particles, remove water.
- Example of saving results (Saving the Results).

## Prepare the working folder

Before running the scripts, create the root folder where all results will be saved (regardless of where you run the jobs) with <a href="{{ site.baseurl }}/episodes/scripts/course_paths.py">course_paths.py</a>:

```bash
python docs/episodes/scripts/course_paths.py
```

This creates the main folder (default `~/Concepcion26`), creates `data/` and `results/`, and shows the `export COURSE_DIR="..."` you must run in the terminal so all scripts use the same destination.

Then copy the data into `data/`:

```bash
cp -R Course-MD-Data-main/* "$COURSE_DIR/data"
```

## Alanine dipeptide

### Guided demo

<!-- sync-from: docs/episodes/scripts/02-preparacion-sistema_simple.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/02-preparacion-sistema_simple.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/02-preparacion-sistema_simple.ipynb" download>Download notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/02-preparacion-sistema_simple.py" download>Download script (.py)</a></div></div>

### Exercise

- Run the simple script and review `"$COURSE_DIR/results/02-preparacion-sistema/simple/alanine_solvated.pdb"`.

### Notebooks and scripts

- This notebook prepares the alanine system with Modeller, adds solvent, and writes the simple solvated PDB ready for simulation. (<a href="{{ site.baseurl }}/episodes/notebooks/02-preparacion-sistema_simple.ipynb">notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/02-preparacion-sistema_simple.py">script</a>)

## Protein-ligand complex

### Guided demo

<!-- sync-from: docs/episodes/scripts/02-preparacion-sistema.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/02-preparacion-sistema.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/02-preparacion-sistema.ipynb" download>Download notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/02-preparacion-sistema.py" download>Download script (.py)</a></div></div>

### Exercise

- Run the complex script with `--output-dir outputs`.
- Try a different `--output-base`.

### Notebooks and scripts

- This notebook walks through Modeller/PDBFixer to repair the protein-ligand complex, allowing you to add solvent, ions, and save preparatory files. (<a href="{{ site.baseurl }}/episodes/notebooks/02-preparacion-sistema.ipynb">notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/02-preparacion-sistema.py">script</a>)

## Key points

- Preparation depends on the target system.
- Saving prepared versions improves traceability.

## Modeller (OpenMM guide)

### Add hydrogens

<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/02-preparacion-sistema.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/02-preparacion-sistema.ipynb" download>Download notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/02-preparacion-sistema.py" download>Download script (.py)</a></div></div>

### Add solvent

<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/02-preparacion-sistema.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/02-preparacion-sistema.ipynb" download>Download notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/02-preparacion-sistema.py" download>Download script (.py)</a></div></div>

### Add membrane

<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/02-preparacion-sistema.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/02-preparacion-sistema.ipynb" download>Download notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/02-preparacion-sistema.py" download>Download script (.py)</a></div></div>

### Add extra particles

<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/02-preparacion-sistema.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/02-preparacion-sistema.ipynb" download>Download notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/02-preparacion-sistema.py" download>Download script (.py)</a></div></div>

### Remove water

<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/02-preparacion-sistema.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/02-preparacion-sistema.ipynb" download>Download notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/02-preparacion-sistema.py" download>Download script (.py)</a></div></div>

### Save results (Saving the Results)

<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/02-preparacion-sistema.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/02-preparacion-sistema.ipynb" download>Download notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/02-preparacion-sistema.py" download>Download script (.py)</a></div></div>

Course script: <a href="{{ site.baseurl }}/episodes/scripts/openmm_modeller_save.py">openmm_modeller_save.py</a>

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/03-simulaciones-clasicas/">Previous</a>
  <a href="{{ site.baseurl }}/episodes/">All episodes</a>
  <a href="{{ site.baseurl }}/episodes/05-muestreo-avanzado/">Next</a>
</div>
