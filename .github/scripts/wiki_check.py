#!/usr/bin/env python3
"""Validate the tracked Markdown source for the FusionPep GitHub Wiki.

Adapted from the LangeLab PXAudit wiki check. FusionPep's wiki also carries
images, so an ``images/`` directory is allowed and image links must resolve.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from urllib.parse import unquote

_REQUIRED_PAGES = ("Home", "_Sidebar", "_Footer")
_IMAGE_DIR = "images"
_IMAGE_SUFFIXES = (".png", ".svg", ".jpg", ".jpeg", ".gif", ".webp")
_PORTABLE_FILENAME_CHARS = re.compile(r'[\\/:*?"<>|]')
_WIKI_LINK = re.compile(r"\[\[([^]\n]+)\]\]")
_MARKDOWN_LINK = re.compile(r"(!?)\[([^]\n]*)\]\(([^)\n]+)\)")
_HTML_IMAGE = re.compile(r"""<img\b[^>]*\bsrc\s*=\s*["']([^"']+)["']""", re.IGNORECASE)
_LINK_TITLE = re.compile(r"""\s+(?:"[^"]*"|'[^']*'|\([^)]*\))\s*$""")
_INLINE_CODE = re.compile(r"`+[^`\n]*`+")
_FENCE = re.compile(r"^( {0,3})(`{3,}|~{3,})")


def _blank(text: str) -> str:
    """Replace non-newline characters with spaces so match positions survive."""
    return "".join("\n" if char == "\n" else " " for char in text)


def _without_code(text: str) -> str:
    """Mask fenced and inline code, where link-looking text is not a link."""
    masked: list[str] = []
    fence_char: str | None = None
    fence_length = 0
    for line in text.splitlines(keepends=True):
        match = _FENCE.match(line)
        if fence_char is None:
            if match is not None:
                fence_char = match.group(2)[0]
                fence_length = len(match.group(2))
                masked.append(_blank(line))
            else:
                masked.append(line)
            continue
        masked.append(_blank(line))
        if re.match(rf"^{re.escape(fence_char)}{{{fence_length},}}[ \t]*\r?$", line.lstrip(" ")):
            fence_char = None
            fence_length = 0
    return _INLINE_CODE.sub(lambda match: _blank(match.group(0)), "".join(masked))


def _page_slug(value: str) -> str:
    """Return the GitHub Wiki-equivalent key for a page name."""
    return re.sub(r"[\s_-]+", "-", unquote(value).strip().casefold())


def _line_number(text: str, position: int) -> int:
    return text.count("\n", 0, position) + 1


def _clean_target(raw_target: str) -> str:
    """Strip angle brackets, a link title, a fragment, and a leading ./ from a target."""
    target = unquote(raw_target.strip())
    if target.startswith("<") and ">" in target:
        target = target[1 : target.index(">")]
    else:
        target = _LINK_TITLE.sub("", target).strip()
    return target.split("#", 1)[0].strip().removeprefix("./")


def _is_external(target: str) -> bool:
    return bool(re.match(r"^(?:[a-z][a-z0-9+.-]*:|//)", target, flags=re.IGNORECASE)) or target.startswith("#")


def _check_page_target(target: str, aliases: dict[str, Path], where: str, errors: list[str]) -> None:
    if not target:
        return
    if target.casefold().endswith(".md"):
        errors.append(f"{where} uses a repository-style .md wiki link: {target}")
        target = target[:-3]
    if "/" in target:
        errors.append(f"{where} links outside the wiki: {target}")
    elif _page_slug(target) not in aliases:
        errors.append(f"{where} links to a missing wiki page: {target}")


def _check_image_target(target: str, wiki_dir: Path, where: str, errors: list[str]) -> None:
    if not target.startswith(f"{_IMAGE_DIR}/"):
        errors.append(f"{where} image must be stored under {_IMAGE_DIR}/: {target}")
    elif not (wiki_dir / target).is_file():
        errors.append(f"{where} links to a missing image: {target}")


def _check_images(image_dir: Path, errors: list[str]) -> None:
    for entry in sorted(image_dir.iterdir()):
        if entry.is_symlink() or not entry.is_file():
            errors.append(f"wiki images must be regular files: {entry}")
        elif entry.suffix.casefold() not in _IMAGE_SUFFIXES:
            errors.append(f"unsupported wiki image type: {entry}")
        elif entry.stat().st_size == 0:
            errors.append(f"wiki image is empty: {entry}")


def validate(wiki_dir: Path) -> list[str]:
    """Return validation errors for ``wiki_dir``; an empty list means valid."""
    if not wiki_dir.is_dir():
        return [f"wiki directory not found: {wiki_dir}"]

    errors: list[str] = []
    pages: list[Path] = []
    aliases: dict[str, Path] = {}
    for entry in sorted(wiki_dir.iterdir(), key=lambda path: path.name.casefold()):
        if entry.name == _IMAGE_DIR and entry.is_dir() and not entry.is_symlink():
            _check_images(entry, errors)
            continue
        if entry.is_symlink() or not entry.is_file():
            errors.append(f"wiki entries must be Markdown files or the {_IMAGE_DIR}/ directory: {entry}")
            continue
        if entry.suffix != ".md":
            errors.append(f"wiki pages must be Markdown files: {entry}")
            continue
        if _PORTABLE_FILENAME_CHARS.search(entry.name):
            errors.append(f"wiki filename contains a non-portable character: {entry}")
        if entry.stat().st_size == 0:
            errors.append(f"wiki page is empty: {entry}")
        pages.append(entry)
        for alias in {entry.stem.casefold(), _page_slug(entry.stem)}:
            previous = aliases.setdefault(alias, entry)
            if previous != entry:
                errors.append(f"case-insensitive or GitHub Wiki slug collision: {previous} and {entry}")

    for required in _REQUIRED_PAGES:
        if required.casefold() not in aliases:
            errors.append(f"required wiki page is missing: {wiki_dir / (required + '.md')}")

    for page in pages:
        text = page.read_text(encoding="utf-8")
        visible = _without_code(text)
        for match in _WIKI_LINK.finditer(visible):
            target = _clean_target(match.group(1).split("|", 1)[0])
            if target and not _is_external(target):
                _check_page_target(target, aliases, f"{page}:{_line_number(text, match.start())}", errors)
        for match in _MARKDOWN_LINK.finditer(visible):
            is_image, raw_target = match.group(1) == "!", match.group(3)
            if _is_external(raw_target.strip()):
                continue
            target = _clean_target(raw_target)
            where = f"{page}:{_line_number(text, match.start())}"
            if is_image or target.casefold().endswith(_IMAGE_SUFFIXES):
                _check_image_target(target, wiki_dir, where, errors)
            else:
                _check_page_target(target, aliases, where, errors)
        for match in _HTML_IMAGE.finditer(visible):
            if not _is_external(match.group(1)):
                where = f"{page}:{_line_number(text, match.start())}"
                _check_image_target(_clean_target(match.group(1)), wiki_dir, where, errors)
    return errors


def main() -> int:
    """Validate the requested wiki directory and return its process status."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wiki_dir", nargs="?", type=Path, default=Path("wiki"))
    args = parser.parse_args()
    errors = validate(args.wiki_dir)
    for error in errors:
        print(error, file=sys.stderr)
    if errors:
        return 1
    page_count = sum(1 for path in args.wiki_dir.iterdir() if path.suffix == ".md")
    print(f"validated {page_count} wiki pages in {args.wiki_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
