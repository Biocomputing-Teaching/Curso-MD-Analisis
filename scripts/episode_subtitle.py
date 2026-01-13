#!/usr/bin/env python3
import sys

SUBTITLES = {
    '01-getting-started': 'Episode 1: Introduction (environment and data)',
    '02-running-simulations': 'Episode 2: Running simulations',
    '03-model-building': 'Episode 3: Preparing and editing the system',
    '04-advanced-examples': 'Episode 4: Advanced simulations and free energy',
    '05-trajectory-analysis': 'Episode 5: Trajectory analysis',
    '06-pyemma': 'Episode 6: Markov models with PyEMMA',
    '07-deeptime': 'Episode 7: Models and spectra with deeptime',
}

if len(sys.argv) != 2:
    raise SystemExit('usage: episode_subtitle.py <slug>')
slug = sys.argv[1]
try:
    print(SUBTITLES[slug])
except KeyError:
    raise SystemExit(f'Unknown episode slug: {slug}')
