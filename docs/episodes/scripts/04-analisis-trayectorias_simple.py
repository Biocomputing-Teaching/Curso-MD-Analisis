#!/usr/bin/env python3
import argparse
import os
from pathlib import Path

import mdtraj as md
import matplotlib.pyplot as plt


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
