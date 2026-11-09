"""Build a blind verification sheet for an independent reviewer.

Panels are numbered and carry no captions, because a caption naming the compound
turns the test into a reading exercise. The answer key is written separately so
the reviewer can be given the image alone.
"""

from __future__ import annotations

import sys
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from chemdraft.projection import draw_fischer, draw_newman  # noqa: E402

OUT = ROOT / "out" / "verify"
OUT.mkdir(parents=True, exist_ok=True)

PANELS = [
    ("Newman, butane, anti (180\u00b0)",
     lambda: draw_newman("CCCC", bond=[1, 2], dihedral=180)),
    ("Newman, butane, gauche (60\u00b0)",
     lambda: draw_newman("CCCC", bond=[1, 2], dihedral=60)),
    ("Newman, butane, eclipsed (0\u00b0)",
     lambda: draw_newman("CCCC", bond=[1, 2], dihedral=0)),
    ("Fischer, L-alanine (S, NH2 left)",
     lambda: draw_fischer("C[C@H](N)C(=O)O")),
    ("Fischer, D-alanine (R, NH2 right)",
     lambda: draw_fischer("C[C@@H](N)C(=O)O")),
    ("Fischer, D-glucose (2R,3S,4R,5R)",
     lambda: draw_fischer("O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO")),
]

COLS = 3
CELL_W, CELL_H = 460, 470


def main() -> int:
    tiles = []
    for number, (name, make) in enumerate(PANELS, start=1):
        projection = make()
        path = OUT / f"{number:02d}.png"
        import pymupdf

        doc = pymupdf.open(stream=projection.svg.encode(), filetype="svg")
        doc[0].get_pixmap(matrix=pymupdf.Matrix(2, 2), alpha=False).save(str(path))
        tiles.append((number, Image.open(path).convert("RGB")))
        print(f"  {number}. {name}")

    rows = (len(tiles) + COLS - 1) // COLS
    sheet = Image.new("RGB", (CELL_W * COLS, CELL_H * rows), "white")
    draw = ImageDraw.Draw(sheet)
    try:
        font = ImageFont.truetype("segoeuisb.ttf", 22)
    except OSError:
        font = ImageFont.load_default()

    for index, (number, image) in enumerate(tiles):
        x, y = (index % COLS) * CELL_W, (index // COLS) * CELL_H
        image.thumbnail((CELL_W - 50, CELL_H - 80))
        sheet.paste(image, (x + (CELL_W - image.width) // 2, y + 40))
        draw.text((x + 16, y + 12), str(number), fill="#888888", font=font)

    sheet_path = OUT / "sheet.png"
    sheet.save(sheet_path)

    key = ["# Answer key", "", "Do not send this with the image.", "",
           "| # | expected reading |", "|---|---|"]
    for number, (name, _) in enumerate(PANELS, start=1):
        key.append(f"| {number} | {name} |")
    (OUT / "answer-key.md").write_text("\n".join(key) + "\n", encoding="utf-8")

    print(f"\nsheet  -> {sheet_path}")
    print(f"key    -> {OUT / 'answer-key.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
