---
layout: default
title: Setup
permalink: /setup/
---

## Quick install (OpenMM)

OpenMM can be installed with `conda` or `pip`. The [official guide](https://docs.openmm.org/latest/userguide/application/01_getting_started.html) recommends Miniconda for an isolated, reproducible environment.

### Recommended option (conda)

```bash
conda create -n md-openmm python=3.10
conda activate md-openmm
conda install -c conda-forge openmm
```

If you have an NVIDIA GPU and want CUDA:

```bash
conda install -c conda-forge openmm cuda-version=12
```

### Alternative (pip)

```bash
pip install openmm
```

If you use an NVIDIA GPU with CUDA 12:

```bash
pip install "openmm[cuda12]"
```

If you use an AMD GPU (HIP):

```bash
pip install "openmm[hip6]"
```

### Quick check

```bash
python - <<'PY'
import openmm
print('openmm', openmm.__version__)
PY
```
## OpenMM-Setup (GUI)

The `openmm-setup` app generates ready-to-run scripts and fixes common structural issues.

```bash
conda install -c conda-forge openmm-setup
```

To launch it:

```bash
openmm-setup
```

## Course packages

Inside the `md-openmm` environment, install the packages we use in the episodes:

```bash
conda install -c conda-forge jupyterlab mdanalysis mdtraj deeptime openff-toolkit openmmforcefields pdbfixer
```

PyEMMA is optional (the project is frozen). If you need it:

```bash
pip install pyemma
```



## Course data

We use an external directory controlled by `COURSE_DIR`:

```bash
export COURSE_DIR=~/Concepcion26
if [ ! -d "$COURSE_DIR" ]; then
  mkdir -p "$COURSE_DIR" 
fi
```

Download data:

```bash
cd "$COURSE_DIR"
curl -L -o Course-MD-Data.zip https://github.com/Biocomputing-Teaching/Course-MD-Data/archive/refs/heads/main.zip
unzip -q Course-MD-Data.zip
rm -rf "$COURSE_DIR/data"
mkdir -p "$COURSE_DIR/data"
cp -R Course-MD-Data-main/* "$COURSE_DIR/data"
rm -rf Course-MD-Data-main Course-MD-Data.zip
```

Download scripts and notebooks:

```bash
cd "$COURSE_DIR"
curl -L -o Curso-MD-Analisis.zip https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/archive/refs/heads/main.zip
unzip -q Curso-MD-Analisis.zip
mkdir -p "$COURSE_DIR/scripts" "$COURSE_DIR/notebooks"
find Curso-MD-Analisis-main/docs/episodes/scripts -type f -name "*.py" -exec cp -i {} "$COURSE_DIR/scripts" \;
find Curso-MD-Analisis-main/docs/episodes/notebooks -type f -name "*.ipynb" -exec cp -i {} "$COURSE_DIR/notebooks" \;
rm -rf Curso-MD-Analisis-main Curso-MD-Analisis.zip
```

Expected structure:

- `data/` input files
- `results/` outputs and analysis
- `scripts/` course scripts
- `notebooks/` course notebooks

## Run simulations (OpenMM scripts)

You can also download the official OpenMM repository to access the example scripts:

```bash
cd "$COURSE_DIR"
git clone git@github.com:openmm/openmm.git
cd openmm/examples/python-examples
```

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/">All episodes</a>
  <a href="{{ site.baseurl }}/episodes/01-introduccion/">Next</a>
</div>
