# /// script
# requires-python = ">=3.12"
# dependencies = ["playwright==1.63.0", "pillow==12.3.0"]
# ///
"""Capture report screenshots and web-sized figures for the README and wiki.

Usage: uv run .github/scripts/demo_images.py OUTPUT_DIR IMAGE_DIR

OUTPUT_DIR is a FusionPep run of the bundled example. The images written to
IMAGE_DIR are the ones the README and wiki reference; build_demo.sh runs this
script after generating the example. Chromium must be installed for Playwright
(`uv run --with playwright==1.63.0 playwright install chromium`).
"""

from __future__ import annotations

import sys
from pathlib import Path

from PIL import Image
from playwright.sync_api import Page, sync_playwright

DESKTOP = {"width": 1440, "height": 900}
PHONE = {"width": 390, "height": 844}
MARGIN = 24  # White space kept around each section capture.
SECTION_MAX_HEIGHT = 1200  # Keep tall sections to a readable crop.
FIGURE_WIDTH = 1600  # Web copy of the 180 mm, 600 dpi figure exports.
FIGURES = ("peptide_coverage", "junction_evidence", "alignment_status")
# Report section id -> collapsible panels to open before the capture.
SECTIONS = {
    "junction": [],
    "peptides": [],
    "coverage": [],
    "sequences": ["Fusion ("],
    "alignment": [],
    "audit": ["Output manifest"],
}


def save_png(image: Image.Image, path: Path) -> None:
    image.save(path, format="PNG", optimize=True)
    print(f"wrote {path} ({image.width}x{image.height}, {path.stat().st_size // 1024} KiB)")


def open_panels(page: Page, section: str, summaries: list[str]) -> None:
    for text in summaries:
        page.locator(f"section#{section} details > summary", has_text=text).first.click()


def capture_section(page: Page, section: str, path: Path) -> None:
    element = page.locator(f"section#{section}")
    element.scroll_into_view_if_needed()
    box = element.bounding_box()
    if box is None:
        raise SystemExit(f"report section #{section} is not visible")
    top = box["y"] + page.evaluate("window.scrollY") - MARGIN
    page.screenshot(
        path=str(path),
        full_page=True,
        clip={
            "x": box["x"] - MARGIN,
            "y": top,
            "width": box["width"] + 2 * MARGIN,
            "height": min(box["height"], SECTION_MAX_HEIGHT) + 2 * MARGIN,
        },
    )
    with Image.open(path) as image:
        save_png(image.copy(), path)


def main() -> int:
    if len(sys.argv) != 3:
        print(__doc__, file=sys.stderr)
        return 2
    output_dir, image_dir = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
    report = output_dir / "fusion_peptide_mapper_report.html"
    if not report.is_file():
        raise SystemExit(f"report not found: {report}")
    image_dir.mkdir(parents=True, exist_ok=True)

    with sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        desktop = browser.new_page(viewport=DESKTOP)
        desktop.goto(report.as_uri(), wait_until="networkidle")
        # The overview runs from the page top to the end of the junction figure:
        # navigation, header, summary cards, and the main decision at a glance.
        figure = desktop.locator("section#junction img").first.bounding_box()
        if figure is None:
            raise SystemExit("the junction figure is not visible in the report")
        overview = image_dir / "report-overview.png"
        desktop.screenshot(
            path=str(overview),
            full_page=True,
            clip={"x": 0, "y": 0, "width": DESKTOP["width"], "height": figure["y"] + figure["height"] + 8},
        )
        with Image.open(overview) as image:
            save_png(image.copy(), overview)
        for section, summaries in SECTIONS.items():
            open_panels(desktop, section, summaries)
            capture_section(desktop, section, image_dir / f"report-{section}.png")

        phone = browser.new_page(viewport=PHONE, device_scale_factor=2, is_mobile=True)
        phone.goto(report.as_uri(), wait_until="networkidle")
        phone_path = image_dir / "report-phone.png"
        phone.screenshot(path=str(phone_path))
        with Image.open(phone_path) as image:
            save_png(image.copy(), phone_path)
        browser.close()

    for name in FIGURES:
        source = output_dir / "figures" / f"{name}.png"
        with Image.open(source) as image:
            height = round(image.height * FIGURE_WIDTH / image.width)
            web = image.convert("RGB").resize((FIGURE_WIDTH, height), Image.Resampling.LANCZOS)
        save_png(web, image_dir / f"figure-{name.replace('_', '-')}.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
