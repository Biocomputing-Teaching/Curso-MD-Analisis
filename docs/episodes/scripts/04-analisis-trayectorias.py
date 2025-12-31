import os
import MDAnalysis as mda
from MDAnalysis.analysis import rms
import matplotlib.pyplot as plt

pdb_path = '../../data/alanine-dipeptide.pdb'
dcd_path = '../../data/alanine-dipeptide.dcd'

if os.path.exists(dcd_path):
    u = mda.Universe(pdb_path, dcd_path)
    output_name = 'rmsd_dcd.png'
else:
    u = mda.Universe('../../data/alanine-dipeptide-multi.pdb')
    output_name = 'rmsd_multimodel.png'

atoms = u.select_atoms('all')
R = rms.RMSD(atoms, atoms)
R.run()

plt.plot(R.results.rmsd[:, 1], R.results.rmsd[:, 2])
plt.xlabel('Frame')
plt.ylabel('RMSD (A)')
plt.title('RMSD')
plt.savefig(output_name, dpi=150)
