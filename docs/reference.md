---
layout: default
title: Referencia
permalink: /reference/
---

## Referencia rápida

### Comandos útiles

```bash
# Activar entorno
conda activate md-openmm

# Abrir JupyterLab
jupyter lab
```

### Estructura sugerida del curso

- `docs/episodes/` episodios Carpentries
- `docs/data/` datos de ejemplo
- `docs/figures/` figuras

## Guía para instructores

- Revisar que OpenMM y dependencias funcionen en el aula.
- Preparar un conjunto mínimo de datos de práctica.
- Validar que los notebooks se ejecuten en GPU y en CPU.

## Materiales

- Notebooks de ejemplo en `docs/episodes/`.
- Figuras en `docs/figures/`.

### Notas

- Usa semillas reproducibles cuando el ejercicio lo permita.
- Guarda parámetros de simulación junto con los resultados.

### Verificación de enlaces

```bash
python scripts/check_links.py
```

## Generar trayectoria DCD

Para crear un DCD corto con OpenMM:

```bash
python scripts/generate_example_dcd.py
```

Esto genera `docs/data/alanine-dipeptide.dcd`.
