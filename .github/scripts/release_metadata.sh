#!/usr/bin/env bash
# Check that DESCRIPTION, CITATION.cff, and NEWS.md agree on the version, and
# print the NEWS.md section for a release version.
#
# A release version is X.Y.Z; a development version is X.Y.Z.9000 (R's
# convention) and needs a "# fusionpep (development version)" NEWS section.
#
# Usage: bash release_metadata.sh [release-notes-output]
# Writes version=... and release=true|false to $GITHUB_OUTPUT when set.
set -euo pipefail

notes_path=${1:-/dev/null}
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
  release=false
  [[ "$first_heading" == "# fusionpep (development version)" ]] ||
    fail "development version $version needs '# fusionpep (development version)' as the first NEWS.md heading"
elif [[ "$version" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
  release=true
  [[ "$first_heading" == "# fusionpep $version" ]] ||
    fail "release version $version needs '# fusionpep $version' as the first NEWS.md heading, found '$first_heading'"
  [[ $(grep -c "^# fusionpep $version\$" NEWS.md) -eq 1 ]] || fail "NEWS.md has more than one section for $version"
  grep -Eq '^date-released:[[:space:]]*"?[0-9]{4}-[0-9]{2}-[0-9]{2}"?[[:space:]]*$' CITATION.cff ||
    fail "CITATION.cff needs date-released (YYYY-MM-DD) for release $version"
else
  fail "version '$version' is neither X.Y.Z nor a X.Y.Z.9000 development version"
fi

# The section body runs from the first heading to the next top-level heading.
awk 'NR > 1 && /^# / { exit } NR > 1 { print }' <(sed -n '/^# /,$p' NEWS.md) |
  sed -e '/[^[:space:]]/,$!d' > "$notes_path"
if [[ "$notes_path" != /dev/null && ! -s "$notes_path" ]]; then
  fail "the NEWS.md section for $version is empty"
fi

if [[ -n "${GITHUB_OUTPUT:-}" ]]; then
  {
    echo "version=$version"
    echo "release=$release"
  } >> "$GITHUB_OUTPUT"
fi
echo "version $version (release: $release)"
