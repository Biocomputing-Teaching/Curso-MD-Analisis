from __future__ import annotations

import argparse
import os
from pathlib import Path


DEFAULT_EPISODES = (
    "01-introduccion",
    "02-preparacion-sistema",
    "03-simulaciones-clasicas",
    "04-analisis-trayectorias",
    "05-muestreo-avanzado",
    "06-pyemma",
    "07-deeptime",
)
DEFAULT_SYSTEMS = ("simple", "complex")


def resolve_course_root(override: str | None = None) -> Path:
    base = (override or os.environ.get("COURSE_DIR", "")).strip().strip('"').strip("'")
    if not base:
        base = str(Path.home() / "Concepcion26")
    return Path(os.path.expandvars(base)).expanduser()


def init_course_dirs(root: Path, episodes: tuple[str, ...], systems: tuple[str, ...]) -> None:
    root.mkdir(parents=True, exist_ok=True)
    (root / "data").mkdir(parents=True, exist_ok=True)
    (root / "results").mkdir(parents=True, exist_ok=True)
    for episode in episodes:
        for system in systems:
            (root / "results" / episode / system).mkdir(parents=True, exist_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Initialize the COURSE_DIR structure.")
    parser.add_argument("--course-dir", default=None, help="Override COURSE_DIR for this run.")
    args = parser.parse_args()
    root = resolve_course_root(args.course_dir)
    init_course_dirs(root, DEFAULT_EPISODES, DEFAULT_SYSTEMS)
    print(root)
    print(f'export COURSE_DIR="{root}"')


if __name__ == "__main__":
    main()
