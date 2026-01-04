#!/usr/bin/env python3
import argparse
import os
from pathlib import Path

import mdtraj as md
import plotly.graph_objects as go


def main() -> None:
    parser = argparse.ArgumentParser(description="Analyse complex trajectory")
    parser.add_argument("-p", "--protein", default="output_minimised.pdb", help="Protein PDB file")
    parser.add_argument("-t", "--trajectory", default="output_traj.dcd", help="Trajectory DCD file")
    parser.add_argument("-o", "--output", default="output_reimaged", help="Output base name")
    parser.add_argument("-r", "--remove-waters", action="store_true", help="Remove waters and ions")
    args = parser.parse_args()

    print("Reading trajectory", args.trajectory)
    traj = md.load(args.trajectory, top=args.protein)
    traj.image_molecules(inplace=True)

    if args.remove_waters:
        print("Removing waters")
        traj = traj.atom_slice(traj.top.select("not resname HOH POPC CL NA"))

    print("Realigning")
    prot = traj.top.select("protein")
    traj.superpose(traj[0], atom_indices=prot)

    out_dir = COURSE_DIR / "results" / "04-analisis-trayectorias" / "complex"
    out_dir.mkdir(parents=True, exist_ok=True)
    output_base = str(out_dir / args.output)

    print("Writing re-imaged PDB", output_base + ".pdb")
    traj[0].save(output_base + ".pdb")

    print("Writing re-imaged trajectory", output_base + ".dcd")
    traj.save(output_base + ".dcd")

    atoms = traj.top.select("chainid 1")
    rmsd_lig = md.rmsd(traj, traj, frame=0, atom_indices=atoms, parallel=True, precentered=False)

    atoms = traj.top.select("chainid 0 and backbone")
    rmsd_bck = md.rmsd(traj, traj, frame=0, atom_indices=atoms, parallel=True, precentered=False)

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=traj.time, y=rmsd_lig, mode="lines", name="Ligand"))
    fig.add_trace(go.Scatter(x=traj.time, y=rmsd_bck, mode="lines", name="Backbone"))

    fig.update_layout(title="Trajectory for " + args.trajectory, xaxis_title="Frame", yaxis_title="RMSD")

    file = output_base + ".svg"
    print("Writing RMSD output to", file)
    fig.write_image(file)


if __name__ == "__main__":
    main()
