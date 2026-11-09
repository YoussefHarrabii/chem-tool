"""Scene composition: molecules, arrows, tables, panels, mechanism arrows and
free primitives on one canvas.

Every molecule is drawn and verified by the same code path as a single figure,
then embedded. A scheme cannot contain a structure a standalone render would
have refused.

Layout is declarative. Containers (`row`, `col`, `grid`, `table`, `panel`)
compute boxes; leaf items declare what to draw. Arrows and mechanism arrows
reference items by id, so they follow the layout instead of hard-coded
coordinates.

Primitives (`line`, `polyline`, `polygon`, `circle`, `rect`) take a box size and
draw in local coordinates, which is what makes non-molecular diagrams possible:
energy profiles, axes, orbital boxes, brackets, apparatus.
"""

from __future__ import annotations

import io
import json
import math
import re
from collections import Counter
from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path
from xml.sax.saxutils import escape

import numpy as np
import pymupdf
from PIL import Image, ImageFont
from rdkit.Chem import rdMolDescriptors

from .chemistry import ChemError, parse
from .render import RenderRequest, RenderResult, render

DEFAULT_FONT = 20
DEFAULT_BOND = 46
# Extra viewport beyond the computed layout. Text can render a few percent wider
# than measured, and anything outside the viewport is clipped before we ever get
# a chance to crop to the ink.
SLACK = 48
FONT_STACK = "DejaVu Sans, Segoe UI, Arial, sans-serif"

CONTAINERS = {"row", "col", "grid", "table", "panel", "reaction"}
PRIMITIVES = {"line", "polyline", "polygon", "circle", "rect"}
PROJECTIONS = {"newman", "fischer"}
LEAVES = {"molecule", "text", "arrow", "curved_arrow"} | PRIMITIVES | PROJECTIONS


# ------------------------------------------------------------------ styles ---
@dataclass
class Style:
    font_size: int = DEFAULT_FONT
    bond_px: int = DEFAULT_BOND
    color: bool = False
    reference: bool = True

    @property
    def caption_size(self) -> int:
        return max(12, int(self.font_size * 0.85))

    @property
    def stroke(self) -> float:
        return max(1.2, self.font_size * 0.075)


def _font(size: int) -> ImageFont.FreeTypeFont:
    for name in ("DejaVuSans.ttf", "segoeui.ttf", "arial.ttf"):
        try:
            return ImageFont.truetype(name, size)
        except OSError:
            continue
    return ImageFont.load_default()


def _measure(text: str, size: int) -> tuple[int, int, int]:
    font = _font(size)
    left, top, right, bottom = font.getbbox(text or " ")
    ascent, _descent = font.getmetrics()
    return right - left, bottom - top, ascent


# ------------------------------------------------------------------- items ---
@dataclass
class Box:
    x: int = 0
    y: int = 0
    w: int = 0
    h: int = 0


@dataclass
class Leaf:
    kind: str
    spec: dict
    box: Box = field(default_factory=Box)
    result: RenderResult | None = None
    caption: str = ""
    caption_h: int = 0
    crop_x: int = 0
    crop_y: int = 0
    fragment: str = ""
    notes: list[str] = field(default_factory=list)


@dataclass
class SceneResult:
    svg: str
    png: bytes
    width: int
    height: int
    molecules: int
    caveats: list[str] = field(default_factory=list)


class SpecError(Exception):
    """The scene description cannot be drawn as written."""


# ------------------------------------------------------------------ flatten ---
def _element_counts(smiles_list: list[str]) -> Counter:
    counts: Counter = Counter()
    for smiles in smiles_list:
        formula = rdMolDescriptors.CalcMolFormula(parse(smiles))
        for element, number in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
            counts[element] += int(number or 1)
    return counts


def _check_balance(from_smiles: list[str], to_smiles: list[str],
                   reagents: list[str]) -> None:
    """A reaction arrow asserts that atoms are conserved. Check it, because a
    scheme of individually correct structures can still be one atom short, and
    that kind of error looks perfectly convincing."""
    left = _element_counts(list(from_smiles) + list(reagents))
    right = _element_counts(list(to_smiles))
    if left == right:
        return
    differences = []
    for element in sorted(set(left) | set(right)):
        a, b = left.get(element, 0), right.get(element, 0)
        if a != b:
            differences.append(f"{element}: {a} on the left, {b} on the right")
    raise SpecError(
        "reaction does not balance by element count, so it cannot be drawn: "
        + "; ".join(differences)
        + ". Add the missing reagent or correct a structure."
    )


def _reaction_row(node: dict) -> dict:
    mode = node.get("mode", "skeletal")
    items: list[dict] = []
    for index, smiles in enumerate(node["_from"]):
        if index:
            items.append({"type": "text", "content": "+"})
        items.append({"type": "molecule", "smiles": smiles, "mode": mode})
    items.append({
        "type": "arrow",
        "length": int(node.get("length", 120)),
        "above": str(node.get("above", "")),
        "below": str(node.get("below", "")),
    })
    for index, smiles in enumerate(node["_to"]):
        if index:
            items.append({"type": "text", "content": "+"})
        items.append({"type": "molecule", "smiles": smiles, "mode": mode})
    return {"type": "row", "gap": int(node.get("gap", 14)),
            "align": "center", "items": items}


def _flatten(node: dict, style: Style, leaves: list[Leaf],
             sizes: dict[int, tuple[int, int]],
             plans: dict[int, list[Box]]) -> tuple[int, int]:
    kind = node.get("type")
    if kind not in CONTAINERS and kind not in LEAVES:
        raise SpecError(
            f"unknown item type {kind!r}; use molecule, text, arrow, curved_arrow, "
            "line, polyline, polygon, circle, rect, row, col, grid, table or panel"
        )

    if kind in LEAVES:
        leaf = _build_leaf(node, style)
        leaves.append(leaf)
        sizes[id(node)] = (leaf.box.w, leaf.box.h)
        return leaf.box.w, leaf.box.h

    gap = int(node.get("gap", 24))

    if kind == "reaction":
        from_smiles = node.get("from") or []
        to_smiles = node.get("to") or []
        if not from_smiles or not to_smiles:
            raise SpecError("reaction needs 'from' and 'to' lists of SMILES")
        reagents = node.get("reagents") or []
        _check_balance(list(from_smiles), list(to_smiles), list(reagents))
        node["_from"] = list(from_smiles)
        node["_to"] = list(to_smiles)
        row = _reaction_row(node)
        node["_row"] = row
        size = _flatten(row, style, leaves, sizes, plans)
        sizes[id(node)] = size
        return size

    if kind == "panel":
        children = node.get("items") or []
        if not children:
            raise SpecError("panel has no items")
        inner = {"type": node.get("layout", "col"), "gap": gap, "items": children}
        inner_size = _flatten(inner, style, leaves, sizes, plans)
        pad = int(node.get("padding", 12))
        title = str(node.get("title", ""))
        title_h = _measure(title, style.caption_size)[1] + 10 if title else 0
        title_w = _text_width(title, style.caption_size) if title else 0
        size = (max(inner_size[0], title_w) + 2 * pad,
                inner_size[1] + 2 * pad + title_h)
        sizes[id(node)] = size
        sizes[id(inner)] = inner_size
        node["_inner"] = inner
        return size

    if kind == "table":
        rows = node.get("rows")
        if not rows:
            raise SpecError("table needs a 'rows' list of lists of items")
        cell_sizes = [[_flatten(cell, style, leaves, sizes, plans) for cell in row]
                      for row in rows]
        cols = max(len(row) for row in rows)
        col_w = [0] * cols
        row_h = [0] * len(rows)
        for r, row in enumerate(cell_sizes):
            for c, (w, h) in enumerate(row):
                col_w[c] = max(col_w[c], w)
                row_h[r] = max(row_h[r], h)
        boxes: list[Box] = []
        y = 0
        for r, row in enumerate(cell_sizes):
            x = 0
            for c, (w, h) in enumerate(row):
                boxes.append(Box(x, y, w, h))
                x += col_w[c] + gap
            y += row_h[r] + gap
        plans[id(node)] = boxes
        size = (sum(col_w) + gap * (cols - 1),
                sum(row_h) + gap * (len(rows) - 1))
        sizes[id(node)] = size
        node["_cells"] = [cell for row in rows for cell in row]
        node["_col_w"] = col_w
        node["_row_h"] = row_h
        return size

    children = node.get("items") or []
    if not children:
        raise SpecError(f"{kind} item has no items")
    child_sizes = [_flatten(c, style, leaves, sizes, plans) for c in children]

    if kind == "row":
        size = (sum(w for w, _ in child_sizes) + gap * (len(children) - 1),
                max(h for _, h in child_sizes))
    elif kind == "col":
        size = (max(w for w, _ in child_sizes),
                sum(h for _, h in child_sizes) + gap * (len(children) - 1))
    else:
        cols = int(node.get("cols", 2))
        if cols < 1:
            raise SpecError("grid needs at least one column")
        rows_n = math.ceil(len(children) / cols)
        cell_w = max(w for w, _ in child_sizes)
        cell_h = max(h for _, h in child_sizes)
        size = (cell_w * cols + gap * (cols - 1),
                cell_h * rows_n + gap * (rows_n - 1))

    sizes[id(node)] = size
    return size


def _build_leaf(node: dict, style: Style) -> Leaf:
    kind = node.get("type")
    if kind == "molecule":
        return _molecule_leaf(node, style)
    if kind == "text":
        return _text_leaf(node, style)
    if kind == "arrow":
        return _arrow_leaf(node, style)
    if kind == "curved_arrow":
        leaf = Leaf(kind="curved_arrow", spec=node)
        leaf.box.w = leaf.box.h = 0
        return leaf
    if kind in PROJECTIONS:
        return _projection_leaf(node, style)
    return _primitive_leaf(node, style)


def _projection_leaf(node: dict, style: Style) -> Leaf:
    from .projection import draw_fischer, draw_newman

    smiles = node.get("smiles")
    if not smiles:
        raise SpecError(f"{node.get('type')} item needs a 'smiles' string")
    bond_px = int(node.get("bond_px", style.bond_px))
    font_px = int(node.get("font_size", style.font_size))
    try:
        if node.get("type") == "newman":
            projection = draw_newman(
                smiles,
                bond=node.get("bond"),
                dihedral=float(node.get("dihedral", 60)),
                bond_px=bond_px,
                font_px=font_px,
            )
        else:
            projection = draw_fischer(smiles, bond_px=bond_px, font_px=font_px)
    except ChemError as exc:
        raise SpecError(f"{node.get('type')} {smiles!r}: {exc}") from exc

    leaf = Leaf(kind=node["type"], spec=node, fragment=projection.svg,
                notes=list(projection.notes))
    leaf.box.w = projection.width
    leaf.box.h = projection.height
    return leaf


def _molecule_leaf(node: dict, style: Style) -> Leaf:
    smiles = node.get("smiles")
    if not smiles:
        raise SpecError("molecule item needs a 'smiles' string")
    request = RenderRequest(
        smiles=smiles,
        mode=node.get("mode", "skeletal"),
        highlight=node.get("highlight"),
        bond_px=int(node.get("bond_px", style.bond_px)),
        color=style.color,
        reference=style.reference,
        rasterize=False,
    )
    try:
        result = render(request)
    except ChemError as exc:
        raise SpecError(f"molecule {smiles!r}: {exc}") from exc
    if not result.report.ok:
        raise SpecError(
            f"molecule {smiles!r} failed verification: "
            + "; ".join(result.report.problems)
        )
    caption = node.get("caption", "")
    caption_h = _measure(caption, style.caption_size)[1] + 8 if caption else 0
    leaf = Leaf(kind="molecule", spec=node, result=result,
                caption=caption, caption_h=caption_h)

    bounds = _svg_bbox(result.svg)
    if bounds is None:
        x0, y0, x1, y1 = 0, 0, result.width, result.height
    else:
        pad = 4
        x0 = max(0, int(bounds[0]) - pad)
        y0 = max(0, int(bounds[1]) - pad)
        x1 = int(bounds[2]) + pad
        y1 = int(bounds[3]) + pad

    caption_w = _text_width(caption, style.caption_size)
    leaf.crop_x, leaf.crop_y = x0, y0
    leaf.box.w = max(1, x1 - x0, caption_w)
    leaf.box.h = max(1, y1 - y0) + caption_h
    return leaf


def _svg_bbox(svg: str) -> tuple[float, float, float, float] | None:
    xs: list[float] = []
    ys: list[float] = []
    for match in re.finditer(r"d=['\"]([^'\"]+)['\"]", svg):
        numbers = [float(v) for v in re.findall(r"-?\d+\.?\d*", match.group(1))]
        xs.extend(numbers[0::2])
        ys.extend(numbers[1::2])
    if not xs or not ys:
        return None
    return min(xs), min(ys), max(xs), max(ys)


def _text_width(text: str, size: int) -> int:
    """Measured width plus a margin. The SVG renders with a font stack that may
    resolve differently from the font used to measure, and a few pixels short
    lets a long label run into the next table column."""
    if not text:
        return 0
    return int(_measure(text, size)[0] * 1.10) + 10


def _text_leaf(node: dict, style: Style) -> Leaf:
    content = str(node.get("content", ""))
    size = int(node.get("size", style.font_size))
    w, h, _ = _measure(content, size)
    leaf = Leaf(kind="text", spec=node)
    leaf.box.w = max(w + 4, _text_width(content, size))
    leaf.box.h = h + 4
    return leaf


def _arrow_leaf(node: dict, style: Style) -> Leaf:
    length = int(node.get("length", max(110, style.font_size * 5)))
    above = str(node.get("above", ""))
    below = str(node.get("below", ""))
    label_h = 0
    if above:
        label_h += _measure(above, style.caption_size)[1]
    if below:
        label_h += _measure(below, style.caption_size)[1]
    leaf = Leaf(kind="arrow", spec=node)
    leaf.box.w = length
    leaf.box.h = max(int(style.font_size * 2.0), label_h + int(style.font_size * 1.6))
    return leaf


def _primitive_leaf(node: dict, style: Style) -> Leaf:
    if "w" not in node or "h" not in node:
        raise SpecError(f"{node.get('type')} item needs 'w' and 'h'")
    leaf = Leaf(kind=node["type"], spec=node)
    leaf.box.w = int(node["w"])
    leaf.box.h = int(node["h"])
    if leaf.box.w <= 0 or leaf.box.h <= 0:
        raise SpecError(f"{node.get('type')} item needs positive w and h")
    return leaf


# --------------------------------------------------------------------- place ---
def _place(node: dict, box: Box, style: Style, leaves: list[Leaf], cursor: int,
           sizes: dict[int, tuple[int, int]], plans: dict[int, list[Box]],
           boxes: dict[int, Box]) -> int:
    kind = node.get("type")
    boxes[id(node)] = box

    if kind in LEAVES:
        leaves[cursor].box.x, leaves[cursor].box.y = box.x, box.y
        return cursor + 1

    gap = int(node.get("gap", 24))

    if kind == "reaction":
        return _place(node["_row"], box, style, leaves, cursor, sizes, plans, boxes)

    if kind == "panel":
        inner = node["_inner"]
        pad = int(node.get("padding", 12))
        title = str(node.get("title", ""))
        title_h = _measure(title, style.caption_size)[1] + 10 if title else 0
        return _place(inner, Box(box.x + pad, box.y + pad + title_h,
                                 sizes[id(inner)][0], sizes[id(inner)][1]),
                      style, leaves, cursor, sizes, plans, boxes)

    if kind == "table":
        cells = node["_cells"]
        cell_boxes = plans[id(node)]
        for cell, cell_box in zip(cells, cell_boxes):
            cursor = _place(cell, Box(box.x + cell_box.x, box.y + cell_box.y,
                                      cell_box.w, cell_box.h),
                            style, leaves, cursor, sizes, plans, boxes)
        return cursor

    children = node["items"]
    align = node.get("align", "center")

    if kind == "row":
        x = box.x
        for child in children:
            w, h = sizes[id(child)]
            y = box.y if align == "top" else (
                box.y + box.h - h if align == "bottom" else box.y + (box.h - h) // 2)
            cursor = _place(child, Box(x, y, w, h), style, leaves, cursor,
                            sizes, plans, boxes)
            x += w + gap
        return cursor

    if kind == "col":
        y = box.y
        for child in children:
            w, h = sizes[id(child)]
            x = box.x if align == "left" else (
                box.x + box.w - w if align == "right" else box.x + (box.w - w) // 2)
            cursor = _place(child, Box(x, y, w, h), style, leaves, cursor,
                            sizes, plans, boxes)
            y += h + gap
        return cursor

    cols = int(node.get("cols", 2))
    cell_w = max(sizes[id(c)][0] for c in children)
    cell_h = max(sizes[id(c)][1] for c in children)
    for i, child in enumerate(children):
        w, h = sizes[id(child)]
        cx = box.x + (i % cols) * (cell_w + gap) + (cell_w - w) // 2
        cy = box.y + (i // cols) * (cell_h + gap) + (cell_h - h) // 2
        cursor = _place(child, Box(cx, cy, w, h), style, leaves, cursor,
                        sizes, plans, boxes)
    return cursor


# ---------------------------------------------------------------------- draw ---
def _inner(svg: str) -> str:
    start = svg.find(">", svg.find("<svg")) + 1
    end = svg.rfind("</svg>")
    body = svg[start:end]
    out = []
    for line in body.splitlines():
        if line.lstrip().startswith("<rect") and ("#FFFFFF" in line or "white" in line):
            continue
        out.append(line)
    return "\n".join(out)


def _text_svg(x: float, y: float, text: str, size: int, anchor: str = "start",
              colour: str = "#000000", italic: bool = False) -> str:
    style = ' font-style="italic"' if italic else ""
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT_STACK}" '
            f'font-size="{size}" fill="{colour}" text-anchor="{anchor}"{style}>'
            f'{escape(text)}</text>')


def _arrow_svg(leaf: Leaf, style: Style) -> list[str]:
    b = leaf.box
    parts: list[str] = []
    margin = 6
    y = b.y + b.h // 2
    x1, x2 = b.x + margin, b.x + b.w - margin
    head = max(7, style.font_size * 0.45)
    parts.append(
        f'<line x1="{x1}" y1="{y}" x2="{x2 - head}" y2="{y}" '
        f'stroke="#000000" stroke-width="{style.stroke:.2f}"/>')
    parts.append(
        f'<polygon points="{x2},{y} {x2 - head},{y - head * 0.55} '
        f'{x2 - head},{y + head * 0.55}" fill="#000000"/>')
    above = str(leaf.spec.get("above", ""))
    below = str(leaf.spec.get("below", ""))
    mid = (x1 + x2) // 2
    size = style.caption_size
    if above:
        parts.append(_text_svg(mid, y - int(size * 0.55), above, size, "middle"))
    if below:
        parts.append(_text_svg(mid, y + int(size * 1.25), below, size, "middle"))
    return parts


def _primitive_svg(leaf: Leaf, style: Style) -> list[str]:
    b = leaf.box
    spec = leaf.spec
    stroke = spec.get("stroke", "#000000")
    width = float(spec.get("stroke_width", style.stroke))
    fill = spec.get("fill", "none")
    dash = spec.get("dash")
    dash_attr = f' stroke-dasharray="{dash}"' if dash else ""

    def P(point) -> tuple[float, float]:
        return b.x + float(point[0]), b.y + float(point[1])

    kind = leaf.kind
    if kind == "line":
        x1, y1 = P((spec.get("x1", 0), spec.get("y1", b.h)))
        x2, y2 = P((spec.get("x2", b.w), spec.get("y2", 0)))
        head = max(6, style.font_size * 0.35)
        parts = [f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
                 f'stroke="{stroke}" stroke-width="{width:.2f}"{dash_attr}/>']
        if spec.get("arrow"):
            angle = math.atan2(y2 - y1, x2 - x1)
            bx, by = x2 - head * math.cos(angle), y2 - head * math.sin(angle)
            px, py = -math.sin(angle) * head * 0.5, math.cos(angle) * head * 0.5
            parts.append(f'<polygon points="{x2:.1f},{y2:.1f} {bx + px:.1f},{by + py:.1f} '
                         f'{bx - px:.1f},{by - py:.1f}" fill="{stroke}"/>')
        return parts

    if kind in {"polyline", "polygon"}:
        points = spec.get("points") or []
        if len(points) < 2:
            raise SpecError(f"{kind} item needs at least two 'points'")
        pts = [P(p) for p in points]
        if spec.get("smooth") and len(pts) > 2:
            path = _smooth_path(pts)
        else:
            path = "M " + " L ".join(f"{x:.1f},{y:.1f}" for x, y in pts)
        if kind == "polygon":
            path += " Z"
        return [f'<path d="{path}" fill="{fill}" stroke="{stroke}" '
                f'stroke-width="{width:.2f}"{dash_attr}/>']

    if kind == "circle":
        r = float(spec.get("r", min(b.w, b.h) / 2))
        cx, cy = P((b.w / 2, b.h / 2))
        return [f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" fill="{fill}" '
                f'stroke="{stroke}" stroke-width="{width:.2f}"{dash_attr}/>']

    if kind == "rect":
        x, y = P((spec.get("x", 0), spec.get("y", 0)))
        w = float(spec.get("w_box", b.w))
        h = float(spec.get("h_box", b.h))
        radius = float(spec.get("radius", 0))
        return [f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" '
                f'rx="{radius:.1f}" fill="{fill}" stroke="{stroke}" '
                f'stroke-width="{width:.2f}"{dash_attr}/>']

    raise SpecError(f"unknown primitive {kind!r}")


def _smooth_path(pts: list[tuple[float, float]]) -> str:
    """Catmull-Rom through the points, emitted as cubic BÃ©ziers."""
    out = [f"M {pts[0][0]:.1f},{pts[0][1]:.1f}"]
    for i in range(len(pts) - 1):
        p0 = pts[i - 1] if i > 0 else pts[i]
        p1, p2 = pts[i], pts[i + 1]
        p3 = pts[i + 2] if i + 2 < len(pts) else p2
        c1 = (p1[0] + (p2[0] - p0[0]) / 6, p1[1] + (p2[1] - p0[1]) / 6)
        c2 = (p2[0] - (p3[0] - p1[0]) / 6, p2[1] - (p3[1] - p1[1]) / 6)
        out.append(f"C {c1[0]:.1f},{c1[1]:.1f} {c2[0]:.1f},{c2[1]:.1f} "
                   f"{p2[0]:.1f},{p2[1]:.1f}")
    return " ".join(out)


def _anchor_point(anchor: dict, by_id: dict[str, Leaf]) -> tuple[float, float] | None:
    leaf = by_id.get(str(anchor.get("item", "")))
    if leaf is None or leaf.kind != "molecule" or leaf.result is None:
        return None
    result = leaf.result
    if "atom" in anchor:
        idx = int(anchor["atom"])
        if not 0 <= idx < result.prepared.GetNumAtoms():
            raise SpecError(f"atom {idx} does not exist in item {anchor['item']!r}")
        point = result.drawer.GetDrawCoords(idx)
    elif "bond" in anchor:
        idx = int(anchor["bond"])
        if not 0 <= idx < result.prepared.GetNumBonds():
            raise SpecError(f"bond {idx} does not exist in item {anchor['item']!r}")
        bond = result.prepared.GetBondWithIdx(idx)
        a = result.drawer.GetDrawCoords(bond.GetBeginAtomIdx())
        b = result.drawer.GetDrawCoords(bond.GetEndAtomIdx())
        point = type(a)((a.x + b.x) / 2, (a.y + b.y) / 2)
    else:
        raise SpecError("curved_arrow endpoint needs 'atom' or 'bond'")

    return (leaf.box.x - leaf.crop_x + point.x,
            leaf.box.y - leaf.crop_y + point.y)


def _curved_arrow_svg(leaf: Leaf, style: Style,
                      by_id: dict[str, Leaf]) -> list[str]:
    spec = leaf.spec
    start = _anchor_point(spec.get("from") or {}, by_id)
    end = _anchor_point(spec.get("to") or {}, by_id)
    if start is None or end is None:
        raise SpecError(
            "curved_arrow 'from' and 'to' must name molecule items that exist"
        )
    x1, y1 = start
    x2, y2 = end
    dx, dy = x2 - x1, y2 - y1
    length = math.hypot(dx, dy) or 1.0
    curve = float(spec.get("curve", 0.35))
    nx, ny = -dy / length, dx / length
    mx, my = (x1 + x2) / 2, (y1 + y2) / 2
    cx, cy = mx + nx * curve * length, my + ny * curve * length

    # Keep the ends clear of the atom labels they point at, or the head lands
    # on top of a letter and the arrow reads as noise.
    standoff = float(spec.get("standoff", max(7.0, style.font_size * 0.55)))
    d1 = math.hypot(cx - x1, cy - y1) or 1.0
    d2 = math.hypot(cx - x2, cy - y2) or 1.0
    x1 += (cx - x1) / d1 * standoff
    y1 += (cy - y1) / d1 * standoff
    x2 += (cx - x2) / d2 * standoff
    y2 += (cy - y2) / d2 * standoff

    head = max(7, style.font_size * 0.42)
    tx, ty = x2 - cx, y2 - cy
    tlen = math.hypot(tx, ty) or 1.0
    ux, uy = tx / tlen, ty / tlen
    bx, by = x2 - head * ux, y2 - head * uy
    px, py = -uy * head * 0.45, ux * head * 0.45

    width = float(spec.get("stroke_width", style.stroke))
    parts = [
        f'<path d="M {x1:.1f},{y1:.1f} Q {cx:.1f},{cy:.1f} {x2:.1f},{y2:.1f}" '
        f'fill="none" stroke="#000000" stroke-width="{width:.2f}"/>',
        f'<polygon points="{x2:.1f},{y2:.1f} {bx + px:.1f},{by + py:.1f} '
        f'{bx - px:.1f},{by - py:.1f}" fill="#000000"/>',
    ]
    label = str(spec.get("label", ""))
    if label:
        t = 0.5
        lx = (1 - t) ** 2 * x1 + 2 * (1 - t) * t * cx + t ** 2 * x2
        ly = (1 - t) ** 2 * y1 + 2 * (1 - t) * t * cy + t ** 2 * y2
        parts.append(_text_svg(lx, ly - 6, label, style.caption_size, "middle"))
    return parts


def _chrome(node: dict, boxes: dict[int, Box], style: Style) -> list[str]:
    """Borders and titles that must sit behind the items they contain."""
    parts: list[str] = []
    kind = node.get("type")
    box = boxes.get(id(node))
    if kind == "panel" and box is not None:
        parts.append(
            f'<rect x="{box.x}" y="{box.y}" width="{box.w}" height="{box.h}" '
            f'rx="6" fill="none" stroke="#999999" '
            f'stroke-width="{max(1.0, style.stroke * 0.8):.2f}"/>')
        title = str(node.get("title", ""))
        if title:
            pad = int(node.get("padding", 12))
            ascent = _measure(title, style.caption_size)[2]
            parts.append(_text_svg(box.x + pad, box.y + 6 + ascent, title,
                                   style.caption_size, colour="#555555"))
    elif kind == "table" and box is not None and node.get("border"):
        gap = int(node.get("gap", 24))
        col_w = node.get("_col_w") or []
        row_h = node.get("_row_h") or []
        x = box.x
        for width in col_w[:-1]:
            x += width + gap / 2
            parts.append(f'<line x1="{x:.1f}" y1="{box.y}" x2="{x:.1f}" '
                         f'y2="{box.y + box.h}" stroke="#999999" stroke-width="1"/>')
            x += gap / 2
        y = box.y
        for height in row_h[:-1]:
            y += height + gap / 2
            parts.append(f'<line x1="{box.x}" y1="{y:.1f}" x2="{box.x + box.w}" '
                         f'y2="{y:.1f}" stroke="#999999" stroke-width="1"/>')
            y += gap / 2

    for child in node.get("items") or []:
        parts.extend(_chrome(child, boxes, style))
    for cell in node.get("_cells") or []:
        parts.extend(_chrome(cell, boxes, style))
    if "_inner" in node:
        parts.extend(_chrome(node["_inner"], boxes, style))
    if "_row" in node:
        parts.extend(_chrome(node["_row"], boxes, style))
    return parts


def _draw_leaf(leaf: Leaf, style: Style, by_id: dict[str, Leaf]) -> list[str]:
    if leaf.kind == "molecule":
        assert leaf.result is not None
        parts = [f'<g transform="translate({leaf.box.x - leaf.crop_x},'
                 f'{leaf.box.y - leaf.crop_y})">', _inner(leaf.result.svg), "</g>"]
        if leaf.caption:
            parts.append(_text_svg(leaf.box.x + leaf.box.w // 2,
                                   leaf.box.y + leaf.box.h - 6, leaf.caption,
                                   style.caption_size, "middle", colour="#333333"))
        return parts
    if leaf.kind == "text":
        size = int(leaf.spec.get("size", style.font_size))
        w, h, ascent = _measure(str(leaf.spec.get("content", "")), size)
        anchor = leaf.spec.get("align", "start")
        x = leaf.box.x + (leaf.box.w - w) // 2 if anchor == "center" else leaf.box.x
        italic = bool(leaf.spec.get("italic", False))
        return [_text_svg(x, leaf.box.y + ascent,
                          str(leaf.spec.get("content", "")), size, "start",
                          italic=italic)]
    if leaf.kind == "arrow":
        return _arrow_svg(leaf, style)
    if leaf.kind == "curved_arrow":
        return _curved_arrow_svg(leaf, style, by_id)
    if leaf.kind in PROJECTIONS:
        return [f'<g transform="translate({leaf.box.x},{leaf.box.y})">',
                _inner(leaf.fragment), "</g>"]
    return _primitive_svg(leaf, style)


# --------------------------------------------------------------------- scene ---
def render_scene(spec: dict) -> SceneResult:
    style_spec = spec.get("style") or {}
    style = Style(
        font_size=int(style_spec.get("font_size", DEFAULT_FONT)),
        bond_px=int(style_spec.get("bond_px", DEFAULT_BOND)),
        color=bool(style_spec.get("color", False)),
        reference=bool(style_spec.get("reference", True)),
    )

    items = spec.get("items")
    if not items:
        raise SpecError("scene needs a non-empty 'items' list")

    leaves: list[Leaf] = []
    sizes: dict[int, tuple[int, int]] = {}
    plans: dict[int, list[Box]] = {}
    boxes: dict[int, Box] = {}

    root = {"type": "col", "gap": 26, "items": items}
    width, height = _flatten(root, style, leaves, sizes, plans)
    _place(root, Box(0, 0, width, height), style, leaves, 0, sizes, plans, boxes)

    by_id: dict[str, Leaf] = {}
    for leaf in leaves:
        ident = leaf.spec.get("id")
        if ident:
            if ident in by_id:
                raise SpecError(f"duplicate item id {ident!r}")
            by_id[str(ident)] = leaf

    padding = int((spec.get("canvas") or {}).get("padding", 24))
    parts: list[str] = _chrome(root, boxes, style)
    caveats: list[str] = []
    for leaf in leaves:
        parts.extend(_draw_leaf(leaf, style, by_id))
        if leaf.kind == "molecule" and leaf.result is not None:
            caveats.extend(leaf.result.report.notes)
        caveats.extend(leaf.notes)
    caveats.extend(_overlap_warnings(leaves))

    total_w = width + 2 * padding + 2 * SLACK
    total_h = height + 2 * padding + 2 * SLACK
    background = (spec.get("canvas") or {}).get("background", "#ffffff")
    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{total_w}" height="{total_h}" '
        f'viewBox="0 0 {total_w} {total_h}">'
        f'<rect x="0" y="0" width="{total_w}" height="{total_h}" fill="{background}"/>'
        f'<g transform="translate({padding},{padding})">'
        + "".join(parts)
        + "</g></svg>"
    )

    doc = pymupdf.open(stream=svg.encode("utf-8"), filetype="svg")
    png = doc[0].get_pixmap(matrix=pymupdf.Matrix(2.0, 2.0), alpha=False).tobytes("png")

    return SceneResult(
        svg=svg, png=_crop(png), width=total_w, height=total_h,
        molecules=sum(1 for leaf in leaves if leaf.kind == "molecule"),
        caveats=sorted(set(caveats)),
    )


def _overlap_warnings(leaves: list[Leaf], margin: int = 4) -> list[str]:
    """Items that share space. Overlapping boxes mean a caption on a bond or two
    structures touching, which reads as one object."""
    warnings: list[str] = []
    for a, b in combinations(leaves, 2):
        if a.box.w == 0 or a.box.h == 0 or b.box.w == 0 or b.box.h == 0:
            continue
        overlap_x = min(a.box.x + a.box.w, b.box.x + b.box.w) - max(a.box.x, b.box.x)
        overlap_y = min(a.box.y + a.box.h, b.box.y + b.box.h) - max(a.box.y, b.box.y)
        if overlap_x > margin and overlap_y > margin:
            warnings.append(
                f"items overlap by {overlap_x}x{overlap_y} px "
                f"({a.kind} and {b.kind}), so something may be drawn on top of "
                f"something else"
            )
    return warnings


def _crop(png: bytes, pad: int = 16) -> bytes:
    with Image.open(io.BytesIO(png)) as raw:
        im = raw.convert("RGB")
    grey = np.asarray(im.convert("L"))
    ink = grey < 250
    if not ink.any():
        return png
    ys, xs = np.where(ink)
    x0, y0 = max(0, int(xs.min()) - pad), max(0, int(ys.min()) - pad)
    x1 = min(grey.shape[1] - 1, int(xs.max()) + pad)
    y1 = min(grey.shape[0] - 1, int(ys.max()) + pad)
    im = im.crop((x0, y0, x1 + 1, y1 + 1))
    buf = io.BytesIO()
    im.save(buf, format="PNG", optimize=True)
    return buf.getvalue()


def load_spec(path: str | Path) -> dict:
    try:
        # utf-8-sig, because a file written by a Windows shell or editor often
        # carries a byte-order mark and plain utf-8 decoding refuses it.
        return json.loads(Path(path).read_text(encoding="utf-8-sig"))
    except FileNotFoundError as exc:
        raise SpecError(f"spec file not found: {path}") from exc
    except json.JSONDecodeError as exc:
        raise SpecError(f"spec is not valid JSON: {exc}") from exc


