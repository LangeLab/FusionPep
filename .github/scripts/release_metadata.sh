#!/usr/bin/env bash
# Check the version metadata in DESCRIPTION, CITATION.cff, and CHANGELOG.md.
#
# Every run checks that DESCRIPTION and CITATION.cff carry the same X.Y.Z
# version and that each CHANGELOG.md section heading is either
# "## [Unreleased]" (first, at most once) or "## [X.Y.Z] - YYYY-MM-DD", with no
# version listed twice.
#
# With --release TAG, also check that TAG is vX.Y.Z for the DESCRIPTION
# version, that CHANGELOG.md has a non-empty section for X.Y.Z whose date
# matches date-released in CITATION.cff, and write that section to --notes.
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
[[ "$version" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]] || fail "DESCRIPTION version '$version' is not X.Y.Z"
cff_version=$(sed -n 's/^version:[[:space:]]*"\{0,1\}\([^"]*\)"\{0,1\}[[:space:]]*$/\1/p' CITATION.cff)
[[ "$cff_version" == "$version" ]] || fail "CITATION.cff version '$cff_version' differs from DESCRIPTION version '$version'"

changelog=CHANGELOG.md
[[ -f "$changelog" ]] || fail "$changelog is missing"
section_pattern='^## \[([0-9]+\.[0-9]+\.[0-9]+)\] - ([0-9]{4}-[0-9]{2}-[0-9]{2})$'
seen_versions=" "
section_count=0
while IFS= read -r heading; do
  section_count=$((section_count + 1))
  if [[ "$heading" == "## [Unreleased]" ]]; then
    [[ "$section_count" -eq 1 ]] || fail "'## [Unreleased]' must be the first changelog section"
  elif [[ "$heading" =~ $section_pattern ]]; then
    listed=${BASH_REMATCH[1]}
    [[ "$seen_versions" != *" $listed "* ]] || fail "$changelog lists version $listed more than once"
    seen_versions+="$listed "
  else
    fail "changelog heading must be '## [Unreleased]' or '## [X.Y.Z] - YYYY-MM-DD': $heading"
  fi
done < <(grep '^## ' "$changelog" || true)
[[ "$section_count" -gt 0 ]] || fail "$changelog has no sections"

if [[ -n "$release_tag" ]]; then
  [[ "$release_tag" =~ ^v[0-9]+\.[0-9]+\.[0-9]+$ ]] || fail "release tag must have the form vX.Y.Z: $release_tag"
  [[ "$release_tag" == "v$version" ]] || fail "tag $release_tag does not match DESCRIPTION version $version"

  heading=$(grep -E "^## \[${version//./\\.}\] - " "$changelog" || true)
  [[ -n "$heading" ]] || fail "$changelog has no '## [$version] - YYYY-MM-DD' section; rename [Unreleased] before tagging"
  [[ "$heading" =~ $section_pattern ]]
  changelog_date=${BASH_REMATCH[2]}
  cff_date=$(sed -n 's/^date-released:[[:space:]]*"\{0,1\}\([0-9-]*\)"\{0,1\}[[:space:]]*$/\1/p' CITATION.cff)
  [[ "$cff_date" == "$changelog_date" ]] ||
    fail "CITATION.cff date-released '$cff_date' differs from the changelog date $changelog_date for $version"

  # The section body: from its heading to the next "## " heading.
  awk -v heading="$heading" '
    $0 == heading { found = 1; next }
    found && /^## / { exit }
    found { print }
  ' "$changelog" | sed -e '/[^[:space:]]/,$!d' > "$notes_path"
  if [[ "$notes_path" != /dev/null && ! -s "$notes_path" ]]; then
    fail "the changelog section for $version is empty"
  fi
fi
echo "version $version: DESCRIPTION, CITATION.cff, and $changelog are consistent${release_tag:+; ready to release $release_tag}"
