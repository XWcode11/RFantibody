#!/bin/bash

set -euo pipefail


download_if_missing() {
    local output_name="$1"
    local url="$2"
    local label="$3"

    if [ -s "$output_name" ]; then
        echo "$label already present, skipping."
        return 0
    fi

    echo "Downloading $label..."
    wget -O "$output_name" "$url"
}


# Get the absolute path of the directory where this script is located.
currdir=$(cd "$(dirname "$0")" && pwd)
weights_dir="${currdir}/../weights"

# Ensure the shared weights directory exists, then fill any missing files.
mkdir -p "$weights_dir"
cd "$weights_dir"

download_if_missing "RFdiffusion_Ab.pt" \
    "https://files.ipd.uw.edu/pub/RFantibody/RFdiffusion_Ab.pt" \
    "RFdiffusion weights"
download_if_missing "ProteinMPNN_v48_noise_0.2.pt" \
    "https://files.ipd.uw.edu/pub/RFantibody/ProteinMPNN_v48_noise_0.2.pt" \
    "ProteinMPNN weights"
download_if_missing "RF2_ab.pt" \
    "https://files.ipd.uw.edu/pub/RFantibody/RF2_ab.pt" \
    "RF2 weights"
download_if_missing "RFab_noframework-nosidechains-5-10-23_trainingparamsadded.pt" \
    "https://zenodo.org/records/17488258/files/RFab_noframework-nosidechains-5-10-23_trainingparamsadded.pt?download=1" \
    "RF2 TCR weights"
download_if_missing "antifold_model.pt" \
    "https://opig.stats.ox.ac.uk/data/downloads/AntiFold/models/model.pt" \
    "AntiFold weights"

echo "All weights are available."
