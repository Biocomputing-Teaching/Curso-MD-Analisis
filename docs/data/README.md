Datos de ejemplo para los ejercicios del curso.

- <a href="https://biocomputing-teaching.github.io/Curso-MD-Analisis/data/alanine-dipeptide.pdb">alanine-dipeptide.pdb</a>: estructura mínima para pruebas con OpenMM.
- <a href="https://biocomputing-teaching.github.io/Curso-MD-Analisis/data/alanine-dipeptide-multi.pdb">alanine-dipeptide-multi.pdb</a>: PDB multi-model para análisis de trayectorias sin DCD.
- <a href="https://biocomputing-teaching.github.io/Curso-MD-Analisis/data/example_rmsd.csv">example_rmsd.csv</a>: serie de RMSD de ejemplo para gráficas.

## Generar trayectoria DCD

Para crear un DCD corto con OpenMM:

```bash
python scripts/generate_example_dcd.py
```

Esto genera <a href="https://biocomputing-teaching.github.io/Curso-MD-Analisis/data/alanine-dipeptide.dcd">alanine-dipeptide.dcd</a>.
