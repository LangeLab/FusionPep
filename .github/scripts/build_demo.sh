#!/usr/bin/env bash
# Refresh the report screenshots and figure images used by the README and wiki.
#
# Runs the bundled example into OUTPUT_DIR, then captures the report sections
# and web-sized figures into IMAGE_DIR (default wiki/images). Run from the
# project root after setup_renv.R, with Chromium installed for Playwright.
#
# Usage: bash .github/scripts/build_demo.sh OUTPUT_DIR [IMAGE_DIR]
set -euo pipefail

output_dir=${1:?usage: build_demo.sh OUTPUT_DIR [IMAGE_DIR]}
image_dir=${2:-wiki/images}
[[ -f run_fusion_mapper.R ]] || { echo "run from the FusionPep project root" >&2; exit 1; }
mkdir -p "$output_dir"
output_dir=$(cd "$output_dir" && pwd)

Rscript run_fusion_mapper.R --output="$output_dir"
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
uv run --quiet "$script_dir/demo_images.py" "$output_dir" "$image_dir"
echo "example output: $output_dir"
echo "images: $image_dir"
