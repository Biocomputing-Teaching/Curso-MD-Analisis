#!/usr/bin/env bash
set -euo pipefail

# ----------------------------
# Download MD figures (Wikipedia/Wikimedia)
# ----------------------------

OUTDIR="/home/jordivilla/GitHub/TEACHING/Curso-MD-Analisis/docs/figures"
mkdir -p "$OUTDIR"

# Prefer curl; fallback to wget
fetch() {
  local url="$1"
  local out="$2"
  if command -v curl >/dev/null 2>&1; then
    curl -L --fail --retry 3 --retry-delay 1 -o "$out" "$url"
  elif command -v wget >/dev/null 2>&1; then
    wget -O "$out" "$url"
  else
    echo "Error: need curl or wget."
    exit 1
  fi
}

# Optional: convert SVG to PDF/PNG if you want LaTeX-friendly formats
svg_to_pdf() {
  local svg="$1"
  local pdf="$2"
  if command -v rsvg-convert >/dev/null 2>&1; then
    rsvg-convert -f pdf -o "$pdf" "$svg"
  elif command -v inkscape >/dev/null 2>&1; then
    inkscape "$svg" --export-type=pdf --export-filename="$pdf" >/dev/null 2>&1
  else
    echo "Note: No SVG converter found (rsvg-convert or inkscape). Keeping SVG: $svg"
    return 0
  fi
}

# ----------------------------
# Define figures
# - Use Wikimedia "Special:FilePath/<Filename>" for stable direct downloads
# - If you later want different images, swap filenames/URLs here.
# ----------------------------

# Each entry: output_basename | direct_download_url | source_page_url
FIGS=(
  # Molecular mechanics potential energy function (image used in Wikipedia topic area)
  "mm_potential_energy|https://commons.wikimedia.org/wiki/Special:FilePath/MM_PEF.png|https://commons.wikimedia.org/wiki/File:MM_PEF.png"

  # Verlet integration (use a representative diagram; may be SVG on Commons)
  # If this particular file ever changes, swap it for another from the Verlet integration article.
  "verlet_integration|https://commons.wikimedia.org/wiki/Special:FilePath/Areal_velocity_Newton4.svg|https://commons.wikimedia.org/wiki/File:Areal_velocity_Newton4.svg"

  # Replica exchange / parallel tempering schematic
"replica_exchange_md|https://commons.wikimedia.org/wiki/Special:FilePath/Schematic_of_a_replica_exchange_molecular_dynamics_simulation.svg||https://commons.wikimedia.org/wiki/File:Schematic_of_a_replica_exchange_molecular_dynamics_simulation.svg"

  # Autocorrelation illustration (use a standard ACF plot from Commons)
  "autocorrelation|https://commons.wikimedia.org/wiki/Special:FilePath/Autokorrelation.png|https://commons.wikimedia.org/wiki/File:Autokorrelation.png"


  # Hidden Markov model trellis (often used in Wikipedia HMM page)
"hmm_trellis|https://commons.wikimedia.org/wiki/Special:FilePath/Hidden_Markov_model.svg|https://commons.wikimedia.org/wiki/File:Hidden_Markov_model.svg"


  # Transfer operator (schematic; choose a Commons file)
"transfer_operator|https://commons.wikimedia.org/wiki/Special:FilePath/Illustration_of_the_Propagation_of_an_Ensemble_by_a_Flow.png|https://commons.wikimedia.org/wiki/File:Illustration_of_the_Propagation_of_an_Ensemble_by_a_Flow.png"


)

# Print what will be downloaded before starting
echo "Planned figure downloads:"
for entry in "${FIGS[@]}"; do
  IFS="|" read -r base url src <<< "$entry"
  echo "  - $base -> $url"
done
echo ""

# ----------------------------
# Download
# ----------------------------
echo "# Figures sources" > FIGURES_SOURCES.md
echo "" >> FIGURES_SOURCES.md
echo "Downloaded from Wikimedia Commons. Each file is linked to its Commons page for attribution/licensing." >> FIGURES_SOURCES.md
echo "" >> FIGURES_SOURCES.md

for entry in "${FIGS[@]}"; do
  IFS="|" read -r base url src <<< "$entry"

  # Infer extension from URL (simple heuristic)
  ext="${url##*.}"
  out="$OUTDIR/${base}.${ext}"

  download_status="downloaded"
  if [[ -f "$out" ]]; then
    echo "Skipping download (already exists): $out"
    download_status="cached"
  else
    echo "Downloading: $base -> $out"
    if ! fetch "$url" "$out"; then
      echo "Warning: download failed for $url"
      rm -f "$out" >/dev/null 2>&1 || true
      download_status="failed"
    fi
  fi

  # If SVG, also create a PDF alongside it (recommended for LaTeX)
  if [[ "$ext" == "svg" && "$download_status" != "failed" ]]; then
    pdf="$OUTDIR/${base}.pdf"
    if [[ -f "$pdf" ]]; then
      echo "Skipping SVG conversion (PDF exists): $pdf"
    else
      svg_to_pdf "$out" "$pdf"
    fi
  fi

  {
    echo "## ${base}"
    echo "- Direct download: ${url}"
    if [[ "$download_status" == "failed" ]]; then
      echo "- Download result: FAILED (see script output)"
      echo ""
      continue
    fi
    echo "- Source page: ${src}"
    echo ""
  } >> FIGURES_SOURCES.md
done

echo "Done."
echo "Figures in: $OUTDIR/"
echo "Sources file: FIGURES_SOURCES.md"
