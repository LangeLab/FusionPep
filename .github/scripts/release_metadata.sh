#!/usr/bin/env bash
# Check that DESCRIPTION, CITATION.cff, and NEWS.md agree on the version.
#
# Versions are X.Y.Z, or X.Y.Z.9000 during development (R's convention). The
# first NEWS.md heading must be "# fusionpep X.Y.Z" or, for a development
# version, "# fusionpep (development version)".
#
# With --release TAG, also check that TAG is vX.Y.Z for the DESCRIPTION
# version, that CITATION.cff has date-released, and that the NEWS section has
# notes; the notes are written to --notes PATH.
#
# Usage: bash release_metadata.sh [--release TAG] [--notes PATH]
set -euo pipefail

release_tag=""
notes_path=/dev/null
while [[ $# -gt 0 ]]; do
  case "$1" in
    --release) release_tag=${2:?--release needs a tag}; shift 2 ;;
    --notes) notes_path=${2:?--notes needs a path}; shift 2 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

fail() {
  echo "release metadata: $*" >&2
  exit 1
}

version=$(sed -n 's/^Version:[[:space:]]*//p' DESCRIPTION)
[[ -n "$version" ]] || fail "DESCRIPTION has no Version field"
cff_version=$(sed -n 's/^version:[[:space:]]*"\{0,1\}\([^"]*\)"\{0,1\}[[:space:]]*$/\1/p' CITATION.cff)
[[ "$cff_version" == "$version" ]] || fail "CITATION.cff version '$cff_version' differs from DESCRIPTION version '$version'"
[[ -f NEWS.md ]] || fail "NEWS.md is missing"

first_heading=$(grep -m 1 '^# ' NEWS.md || true)
if [[ "$version" =~ ^[0-9]+\.[0-9]+\.[0-9]+\.9[0-9]{3}$ ]]; then
  [[ "$first_heading" == "# fusionpep (development version)" ]] ||
    fail "development version $version needs '# fusionpep (development version)' as the first NEWS.md heading"
elif [[ "$version" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
  [[ "$first_heading" == "# fusionpep $version" ]] ||
    fail "version $version needs '# fusionpep $version' as the first NEWS.md heading, found '$first_heading'"
  [[ $(grep -c "^# fusionpep $version\$" NEWS.md) -eq 1 ]] || fail "NEWS.md has more than one section for $version"
else
  fail "version '$version' is neither X.Y.Z nor a X.Y.Z.9000 development version"
fi

if [[ -n "$release_tag" ]]; then
  [[ "$release_tag" =~ ^v[0-9]+\.[0-9]+\.[0-9]+$ ]] || fail "release tag must have the form vX.Y.Z: $release_tag"
  [[ "$release_tag" == "v$version" ]] || fail "tag $release_tag does not match DESCRIPTION version $version"
  grep -Eq '^date-released:[[:space:]]*"?[0-9]{4}-[0-9]{2}-[0-9]{2}"?[[:space:]]*$' CITATION.cff ||
    fail "releasing $version needs date-released (YYYY-MM-DD) in CITATION.cff"
  # The first NEWS.md section, without its heading or leading blank lines.
  sed -n '/^# /,$p' NEWS.md | awk 'NR > 1 && /^# / { exit } NR > 1 { print }' |
    sed -e '/[^[:space:]]/,$!d' > "$notes_path"
  if [[ "$notes_path" != /dev/null && ! -s "$notes_path" ]]; then
    fail "releasing $version needs notes under '# fusionpep $version' in NEWS.md"
  fi
fi
echo "version $version: DESCRIPTION, CITATION.cff, and NEWS.md agree${release_tag:+; ready to release $release_tag}"
