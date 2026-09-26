#!/usr/bin/env bash
# Build the FusionPep demo from the bundled example.
#
# Runs the example into SITE_DIR (the full output folder, with an index page
# that opens the report) and captures the report screenshots and web-sized
# figures into IMAGE_DIR, which defaults to wiki/images: the README and wiki
# reference those files. Run from the project root after setup_renv.R.
#
# Usage: bash .github/scripts/build_demo.sh SITE_DIR [IMAGE_DIR]
set -euo pipefail

site_dir=${1:?usage: build_demo.sh SITE_DIR [IMAGE_DIR]}
image_dir=${2:-wiki/images}
[[ -f run_fusion_mapper.R ]] || { echo "run from the FusionPep project root" >&2; exit 1; }
mkdir -p "$site_dir"
site_dir=$(cd "$site_dir" && pwd)

Rscript run_fusion_mapper.R --output="$site_dir"
cat > "$site_dir/index.html" <<'HTML'
<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta http-equiv="refresh" content="0; url=fusion_peptide_mapper_report.html">
<title>FusionPep example report</title>
</head>
<body>
<p><a href="fusion_peptide_mapper_report.html">Open the FusionPep example report</a>.</p>
</body>
</html>
HTML

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
uv run --quiet "$script_dir/demo_images.py" "$site_dir" "$image_dir"
echo "demo site: $site_dir"
echo "demo images: $image_dir"
