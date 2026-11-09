"""Build a blind identification sheet.

Figures are rendered without captions and numbered. The answer key is written
separately, so a student can be asked to name each one without seeing the
answers. This measures the thing that matters for exam preparation: whether the
drawing is recognisable, not whether it is internally consistent.

    python tools/quiz_sheet.py
"""

from __future__ import annotations

import sys
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from chemdraft.chemistry import ChemError  # noqa: E402
from chemdraft.render import RenderRequest, render  # noqa: E402

OUT = ROOT / "out" / "quiz"
OUT.mkdir(parents=True, exist_ok=True)

MOLECULES = [
    ("ethanol", "CCO", "skeletal"),
    ("dimethyl ether", "COC", "skeletal"),
    ("propan-2-ol", "CC(C)O", "skeletal"),
    ("acetic acid", "CC(=O)O", "skeletal"),
    ("propanone (acetone)", "CC(C)=O", "skeletal"),
    ("benzene", "c1ccccc1", "skeletal"),
    ("phenol", "Oc1ccccc1", "skeletal"),
    ("aniline", "Nc1ccccc1", "skeletal"),
    ("toluene", "Cc1ccccc1", "skeletal"),
    ("cyclohexane", "C1CCCCC1", "skeletal"),
    ("ethyl acetate", "CCOC(C)=O", "skeletal"),
    ("pyridine", "c1ccncc1", "skeletal"),
    ("trans-but-2-ene", "C/C=C/C", "skeletal"),
    ("cis-but-2-ene", "C/C=C\\C", "skeletal"),
    ("alanine (L)", "C[C@H](N)C(=O)O", "skeletal"),
    ("camphor", "CC12CCC(CC1=O)C2(C)C", "skeletal"),
    ("ethanol, all atoms", "CCO", "lewis"),
    ("urea, all atoms", "NC(=O)N", "lewis"),
]

COLS = 4
CELL_W, CELL_H = 420, 380


def main() -> int:
    tiles = []
    key = []
    for number, (name, smiles, mode) in enumerate(MOLECULES, start=1):
        try:
            result = render(RenderRequest(smiles=smiles, mode=mode, target_width=380,
                                          bond_px=38, reference=True))
        except ChemError as exc:
            print(f"  {number:2d}. {name}: refused ({exc})")
            continue
        if not result.report.ok:
            print(f"  {number:2d}. {name}: failed verification")
            continue
        path = OUT / f"{number:02d}.png"
        path.write_bytes(result.png)
        tiles.append((number, Image.open(path).convert("RGB")))
        key.append((number, name, result.description.formula,
                    result.description.canonical, result.depiction.split(" (")[0]))
        print(f"  {number:2d}. {name:22s} {result.description.formula:10s} "
              f"{result.depiction.split(' (')[0]}")

    rows = (len(tiles) + COLS - 1) // COLS
    sheet = Image.new("RGB", (CELL_W * COLS, CELL_H * rows), "white")
    draw = ImageDraw.Draw(sheet)
    try:
        number_font = ImageFont.truetype("segoeuisb.ttf", 22)
    except OSError:
        number_font = ImageFont.load_default()

    for index, (number, im) in enumerate(tiles):
        x, y = (index % COLS) * CELL_W, (index // COLS) * CELL_H
        im.thumbnail((CELL_W - 40, CELL_H - 60))
        sheet.paste(im, (x + (CELL_W - im.width) // 2, y + 34))
        draw.text((x + 14, y + 10), f"{number}", fill="#888888", font=number_font)

    sheet_path = OUT / "sheet.png"
    sheet.save(sheet_path)

    key_lines = ["# Answer key", "",
                 "Do not show this to the person being tested.", "",
                 "| # | compound | formula | canonical SMILES | depiction |",
                 "|---|---|---|---|---|"]
    for number, name, formula, smiles, depiction in key:
        key_lines.append(f"| {number} | {name} | {formula} | `{smiles}` | {depiction} |")
    (OUT / "answer-key.md").write_text("\n".join(key_lines) + "\n", encoding="utf-8")

    print(f"\n{len(tiles)} figures -> {sheet_path}")
    print(f"answer key -> {OUT / 'answer-key.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
