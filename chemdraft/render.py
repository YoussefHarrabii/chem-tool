"""Draw a molecule, verify it, then emit. Order matters: verification runs
before any file exists on disk."""

from __future__ import annotations

import io
from dataclasses import dataclass, field

import numpy as np
import pymupdf
from PIL import Image
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

from .chemistry import ChemError, Description, describe, needs_label, parse
from .layout import compute_layout, score_layout
from .annotate import dots_for, inject
from .reference import apply_coordinates, polycyclic, record_for
from .verify import Report, verify

rdDepictor.SetPreferCoordGen(True)

# A bond is a fixed number of pixels so every figure in a set shares a scale.
DEFAULT_BOND_PX = 46
TARGET_PNG_WIDTH = 1400
MIN_CANVAS = 260


@dataclass
class RenderRequest:
    smiles: str
    mode: str = "skeletal"          # skeletal | lewis
    highlight: str | None = None    # SMARTS
    bond_px: int = DEFAULT_BOND_PX
    target_width: int = TARGET_PNG_WIDTH
    color: bool = False
    pad: float = 0.42               # padding around the ink, in bond lengths
    rasterize: bool = True          # skip PNG for bulk checks
    reference: bool = False         # use a published depiction when available
    lone_pairs: bool | None = None   # None means: on for lewis, off for skeletal


@dataclass
class RenderResult:
    svg: str
    png: bytes
    description: Description
    report: Report
    width: int = 0
    height: int = 0
    notes: list[str] = field(default_factory=list)
    depiction: str = ""
    # Kept so a test can prove the gate has teeth by breaking them on purpose.
    prepared: Chem.Mol | None = None
    drawer: object | None = None


def _span(mol: Chem.Mol) -> tuple[float, float]:
    conf = mol.GetConformer()
    xs, ys = [], []
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        xs.append(p.x)
        ys.append(p.y)
    return (max(xs) - min(xs) if xs else 0.0, max(ys) - min(ys) if ys else 0.0)


def _canvas(mol: Chem.Mol, bond_px: int, margin_px: int) -> tuple[int, int]:
    """Canvas big enough to hold the molecule plus its labels, with room to
    spare. Cropping after the fact is safer than sizing to the last pixel."""
    span_x, span_y = _span(mol)
    scale = bond_px / 1.5
    width = int(span_x * scale) + 2 * margin_px
    height = int(span_y * scale) + 2 * margin_px
    return max(width, MIN_CANVAS), max(height, MIN_CANVAS)


def _draw(mol: Chem.Mol, width: int, height: int, bond_px: int, color: bool,
          highlight_atoms: list[int], highlight_bonds: list[int]):
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = drawer.drawOptions()
    if not color:
        opts.useBWAtomPalette()
    opts.fixedBondLength = bond_px
    opts.bondLineWidth = max(1.6, bond_px / 22.0)
    opts.fixedFontSize = int(round(bond_px * 0.52))
    opts.minFontSize = 6
    opts.maxFontSize = 400
    opts.padding = 0.02
    opts.additionalAtomLabelPadding = 0.06
    if highlight_atoms or highlight_bonds:
        opts.highlightRadius = 0.45
        opts.highlightColour = (0.45, 0.72, 0.95) if color else (0.62, 0.80, 0.62)
    drawer.DrawMolecule(
        mol,
        highlightAtoms=highlight_atoms or None,
        highlightBonds=highlight_bonds or None,
    )
    drawer.FinishDrawing()
    return drawer, drawer.GetDrawingText()


def _crop(png: bytes, pad_px: int, target_width: int) -> bytes:
    """Trim to the ink and add an even margin, so the figure sits centred with
    the same whitespace as every other figure."""
    with Image.open(io.BytesIO(png)) as raw:
        im = raw.convert("RGB")
    gray = np.asarray(im.convert("L"))
    ink = gray < 250
    if not ink.any():
        return png
    ys, xs = np.where(ink)
    x0, x1 = int(xs.min()), int(xs.max())
    y0, y1 = int(ys.min()), int(ys.max())
    x0 = max(0, x0 - pad_px)
    y0 = max(0, y0 - pad_px)
    x1 = min(gray.shape[1] - 1, x1 + pad_px)
    y1 = min(gray.shape[0] - 1, y1 + pad_px)
    im = im.crop((x0, y0, x1 + 1, y1 + 1))
    if target_width > 0 and im.width != target_width:
        h = max(1, round(im.height * target_width / im.width))
        im = im.resize((target_width, h), Image.LANCZOS)
    buf = io.BytesIO()
    im.save(buf, format="PNG", optimize=True)
    return buf.getvalue()


def _render_png(svg: str, canvas_w: int, margin_px: int, target_width: int,
                pad_px: int) -> bytes:
    """Rasterise the verified SVG at a zoom that lands near the target width
    after cropping, so nothing is upscaled and blurred."""
    usable = max(canvas_w - 2 * margin_px, 1)
    zoom = max(1.0, target_width / usable)
    doc = pymupdf.open(stream=svg.encode("utf-8"), filetype="svg")
    pix = doc[0].get_pixmap(matrix=pymupdf.Matrix(zoom, zoom), alpha=False)
    return _crop(pix.tobytes("png"), int(pad_px * zoom), target_width)


def _touches_border(png: bytes, threshold: int = 250) -> bool:
    with Image.open(io.BytesIO(png)) as im:
        a = np.asarray(im.convert("L"))
    band = 1
    return any((e < threshold).any() for e in
               (a[:band, :], a[-band:, :], a[:, :band], a[:, -band:]))


def render(request: RenderRequest) -> RenderResult:
    mol = parse(request.smiles)
    description = describe(mol, request.smiles)

    work = Chem.Mol(mol)
    if request.mode == "lewis":
        work = Chem.AddHs(work)
    elif request.mode != "skeletal":
        raise ChemError(f"unknown mode {request.mode!r}, use skeletal or lewis")

    depiction = "automatic layout"
    reference_note = ""
    laid_out = None
    layout_name, layout_score = "", (0, 0, 0, 0.0)

    # A published depiction beats an invented one, but only for the plain
    # skeletal view: reference records carry heavy atoms only, so a Lewis
    # structure still needs its own layout.
    if request.reference and request.mode == "skeletal":
        ref, why = record_for(work)
        if ref is not None and apply_coordinates(work, ref):
            laid_out = work
            layout_score = score_layout(work)
            depiction = f"reference depiction ({why})"
        else:
            reference_note = why

    if laid_out is None:
        laid_out, layout_name, layout_score = compute_layout(work)

    prepared = rdMolDraw2D.PrepareMolForDrawing(laid_out)
    if request.mode == "lewis":
        for atom in prepared.GetAtoms():
            if not needs_label(atom):
                atom.SetProp("atomLabel", atom.GetSymbol())

    show_pairs = (request.lone_pairs if request.lone_pairs is not None
                  else request.mode == "lewis")

    hi_atoms, hi_bonds = [], []
    if request.highlight:
        pattern = Chem.MolFromSmarts(request.highlight)
        if pattern is None:
            raise ChemError(f"highlight is not valid SMARTS: {request.highlight!r}")
        for match in mol.GetSubstructMatches(pattern):
            hi_atoms.extend(match)
        allowed = set(hi_atoms)
        hi_bonds = [
            b.GetIdx() for b in mol.GetBonds()
            if b.GetBeginAtomIdx() in allowed and b.GetEndAtomIdx() in allowed
        ]

    notes: list[str] = []
    if depiction.startswith("automatic") and polycyclic(work):
        notes.append(
            "depiction is automatic and this skeleton is polycyclic, which is the "
            "case automatic layout often draws in a non-standard projection; "
            "compare with a reference picture before teaching from it"
        )
    if reference_note:
        notes.append(f"reference depiction unavailable: {reference_note}")
    if depiction.startswith("automatic") and (layout_score[0] or layout_score[1] > 0.05):
        notes.append(
            f"layout {layout_name!r} still has {layout_score[0]} crossing bonds "
            f"and crowding worth {layout_score[1]:.2f}, so some bonds may look tight"
        )
    pad_px = request.pad * request.bond_px
    margin = int(request.bond_px * 1.6)
    last = None

    for _attempt in range(4):
        width, height = _canvas(prepared, request.bond_px, margin)
        drawer, svg = _draw(prepared, width, height, request.bond_px,
                            request.color, hi_atoms, hi_bonds)
        if show_pairs:
            svg = inject(svg, dots_for(prepared, drawer, request.bond_px * 0.52))
        report = verify(
            mol, prepared, drawer, svg, width, height,
            label_all=(request.mode == "lewis"),
            expected=description,
        )
        png = _render_png(svg, width, margin, request.target_width, pad_px)
        if request.rasterize and _touches_border(png):
            last = (width, height, svg, png, report)
            margin = int(margin * 1.7)
            notes.append("grew the canvas because the drawing reached the edge")
            continue
        if not request.rasterize:
            png = b""
        return RenderResult(
            svg=svg, png=png, description=description, report=report,
            width=width, height=height, notes=notes, depiction=depiction,
            prepared=prepared, drawer=drawer,
        )

    width, height, svg, png, report = last
    report.problems.append(
        "the drawing keeps reaching the canvas edge, so part of it may be clipped"
    )
    return RenderResult(
        svg=svg, png=png, description=description, report=report,
        width=width, height=height, notes=notes, depiction=depiction,
        prepared=prepared, drawer=drawer,
    )




