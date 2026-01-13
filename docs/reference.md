---
layout: default
title: Reference
permalink: /reference/
---

## Quick reference

### Useful commands

```bash
# Activate environment
conda activate md-openmm

# Launch JupyterLab
jupyter lab
```

### Suggested course structure

- `docs/episodes/` Carpentries episodes
- `docs/data/` sample data
- `docs/figures/` figures

## Instructor guide

- Check that OpenMM and dependencies work in the classroom.
- Prepare a minimal set of practice data.
- Validate that the notebooks run on GPU and CPU.

## Materials

- Example notebooks in `docs/episodes/`.
- Figures in `docs/figures/`.

### Notes

- Use reproducible seeds when the exercise allows it.
- Save simulation parameters alongside the results.

### Link check

```bash
python scripts/check_links.py
```

## Generate a DCD trajectory

To create a short DCD with OpenMM:

```bash
python scripts/generate_example_dcd.py
```

This generates <a href="https://github.com/Biocomputing-Teaching/Course-MD-Data/raw/main/alanine-dipeptide.dcd">alanine-dipeptide.dcd</a>.

## Repository technical documentation

For maintenance details and technical documentation, see the README on GitHub:

- <a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/blob/main/README.md">README.md on GitHub</a>
