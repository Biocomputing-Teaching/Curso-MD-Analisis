---
layout: default
title: Episode 1 - Introduction (environment and data)
permalink: /episodes/01-introduccion/
---

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/">All episodes</a>
  <a href="{{ site.baseurl }}/episodes/03-simulaciones-clasicas/">Next</a>
</div>

<!-- toc:start -->
## Table of contents
- [Duration](#duration)
- [Objectives](#objectives)
- [Content](#content)
- [Preliminary sequence and ligand checks](#preliminary-sequence-and-ligand-checks)
- [Alanine dipeptide](#alanine-dipeptide)
- [Protein-ligand complex](#protein-ligand-complex)
- [Key points](#key-points)
<!-- toc:end -->

## Duration

- **Session:** 45 min
- **Exercises:** 30 min

## Objectives

- Identify input files for alanine and the protein-ligand complex.
- Verify OpenMM and OpenFF Toolkit work in the environment.
- Load base data for both systems following the first-steps workflow.

## Content

- PDB and MOL/SDF as inputs.
- Reading alanine with OpenMM.
- Reading protein and ligand with OpenFF Toolkit.
- Working path `COURSE_DIR` and environment validation.

## Preliminary sequence and ligand checks

Before we launch the OpenMM demos, review two focused notebooks to explore the inputs.

- Sequence exploration: `<a href="{{ site.baseurl }}/episodes/notebooks/01-introduction-sequence-check.ipynb">01-introduction-sequence-check.ipynb</a>` now focuses on the `protein2-ligand2.pdb` coordinates and extracts the amino-acid sequence for the two chains before the simulation steps.
- Ligand exploration: `<a href="{{ site.baseurl }}/episodes/notebooks/01-introduction-ligand-check.ipynb">01-introduction-ligand-check.ipynb</a>` explores the new `$COURSE_DIR/data/complex/ligand2.sdf` entry with RDKit, reports its properties, and sketches the ligand before combining it with OpenMM.


## Alanine dipeptide

### Guided demo

<!-- sync-from: docs/episodes/scripts/01-introduccion_simple.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/01-introduccion_simple.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/01-introduccion_simple.ipynb" download>Download notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/01-introduccion_simple.py" download>Download script (.py)</a></div></div>

### Exercise

- Run the simple script and confirm atom counts.

### Notebooks and scripts

- Recorrer esta versión guiada del episodio muestra la estructura del directorio `COURSE_DIR`, explica cómo leer un PDB con Biopython y prepara el sistema de alanina para el flujo de trabajo básico. (<a href="{{ site.baseurl }}/episodes/notebooks/01-introduccion_simple.ipynb">notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/01-introduccion_simple.py">script</a>)

## Protein-ligand complex

### Visual sanity checks

- Open `$COURSE_DIR/data/complex/protein2-ligand2.pdb` with the new `<a href="{{ site.baseurl }}/episodes/notebooks/01-introduccion-protein2-nglview.ipynb">complex viewer</a>` to inspect the binding pose with `nglview`.
- Slice the same structure with `bash docs/episodes/scripts/01-split-protein2-ligand2.sh` and check the standalone protein in `<a href="{{ site.baseurl }}/episodes/notebooks/01-introduccion-protein2-protein-nglview.ipynb">protein viewer</a>` as well as the ligand in `<a href="{{ site.baseurl }}/episodes/notebooks/01-introduccion-protein2-ligand-nglview.ipynb">ligand viewer</a>`.
- These quick viewers keep you aligned with the complex before the OpenMM preparation steps and complement the protein/ligand notebooks listed above.

### Splitting the complex

Recreate the canonical `protein2.pdb` / `ligand2.pdb` pair any time the `$COURSE_DIR/data/complex` folder changes by running the helper script below (copy/paste from the repository). It assumes the full complex resides in `$COURSE_DIR/data/complex/protein2-ligand2.pdb`:

```bash
#!/usr/bin/env bash
set -euo pipefail

COURSE_DIR="${COURSE_DIR:-$HOME/Concepcion26}"
INPUT_PATH="$COURSE_DIR/data/complex/protein2-ligand2.pdb"
OUTPUT_DIR="$COURSE_DIR/data/complex"

if [[ ! -f "$INPUT_PATH" ]]; then
  echo "Error: missing $INPUT_PATH" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

python - <<PY
from pathlib import Path

input_path = Path("$INPUT_PATH")
output_dir = Path("$OUTPUT_DIR")
protein_path = output_dir / "protein2.pdb"
ligand_path = output_dir / "ligand2.pdb"

protein_lines = []
ligand_lines = []

with input_path.open("r") as fh:
    for line in fh:
        if line.startswith("HETATM") and line[17:20].strip() == "SUB":
            ligand_lines.append(line)
        else:
            protein_lines.append(line)

protein_path.write_text("".join(protein_lines))
ligand_path.write_text("".join(ligand_lines))

print(f"Saved protein {protein_path}")
print(f"Saved ligand  {ligand_path}")
PY
```

Run this once after updating the data bundle so every episode reads the same split files.

### Guided demo

<!-- sync-from: docs/episodes/scripts/01-introduccion.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/01-introduccion.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/01-introduccion.ipynb" download>Download notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/01-introduccion.py" download>Download script (.py)</a></div></div>

### Exercise

- Run the complex script and swap the ligand for <a href="https://github.com/Biocomputing-Teaching/Course-MD-Data/raw/main/complex/ligand1.sdf">ligand1.sdf</a>.
- Save the output to a text file.

### Notebooks and scripts

- Esta versión completa repasa la preparación del complejo proteína-ligando: secuencia rápida con Biopython, visualización con nglview y escritura de archivos Amber para la simulación. (<a href="{{ site.baseurl }}/episodes/notebooks/01-introduccion.ipynb">notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/01-introduccion.py">script</a>)
- The new `nglview` viewers keep the same split seen above but in interactive form: <a href="{{ site.baseurl }}/episodes/notebooks/01-introduccion-protein2-nglview.ipynb">complex viewer</a>, <a href="{{ site.baseurl }}/episodes/notebooks/01-introduccion-protein2-protein-nglview.ipynb">protein viewer</a>, and <a href="{{ site.baseurl }}/episodes/notebooks/01-introduccion-protein2-ligand-nglview.ipynb">ligand viewer</a>.

## Key points

- The course works with a simple and a complex system in parallel.
- OpenMM and OpenFF Toolkit cover reading PDB and MOL/SDF.

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/">All episodes</a>
  <a href="{{ site.baseurl }}/episodes/03-simulaciones-clasicas/">Next</a>
</div>
