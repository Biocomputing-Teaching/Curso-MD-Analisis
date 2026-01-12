#!/usr/bin/env python3
import sys

SUBTITLES = {
    '01-getting-started': 'Episodio 1: Introducción (entorno y datos)',
    '02-running-simulations': 'Episodio 2: Ejecución de simulaciones',
    '03-model-building': 'Episodio 3: Preparación y edición del sistema',
    '04-advanced-examples': 'Episodio 4: Simulaciones avanzadas y energía libre',
    '05-trajectory-analysis': 'Episodio 5: Análisis de trayectorias',
    '06-pyemma': 'Episodio 6: Modelos de Markov con PyEMMA',
    '07-deeptime': 'Episodio 7: Modelos y espectros con Deeptime',
}

if len(sys.argv) != 2:
    raise SystemExit('usage: episode_subtitle.py <slug>')
slug = sys.argv[1]
try:
    print(SUBTITLES[slug])
except KeyError:
    raise SystemExit(f'Unknown episode slug: {slug}')
