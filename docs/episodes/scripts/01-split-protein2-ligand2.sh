#!/usr/bin/env bash
set -euo pipefail

COURSE_DIR="${COURSE_DIR:-$HOME/Concepcion26}"
INPUT_PATH="${COURSE_DIR}/data/complex/protein2-ligand2.pdb"
OUTPUT_DIR="${COURSE_DIR}/data/complex"

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

with input_path.open(encoding="utf-8") as fh:
    for line in fh:
        if line.startswith("HETATM") and line[17:20].strip() == "SUB":
            ligand_lines.append(line)
        else:
            protein_lines.append(line)

protein_path.write_text("".join(protein_lines), encoding="utf-8")
ligand_path.write_text("".join(ligand_lines), encoding="utf-8")

print(f"Saved protein {protein_path} ({len(protein_lines)} lines)")
print(f"Saved ligand  {ligand_path} ({len(ligand_lines)} lines)")
PY
