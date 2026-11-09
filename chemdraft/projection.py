"""Newman and Fischer projections.

Both are conventions that encode stereochemistry in the geometry, so they are
computed rather than templated:

- Newman is a true projection of a real 3D conformer. The molecule is embedded,
  the requested dihedral is set, and the substituents are projected onto the
  plane perpendicular to the view axis.
- Fischer places the main chain vertically and horizontal substituents toward
  the viewer, then *verifies* the result: the arrangement is turned back into
  3D, RDKit re-reads the configuration, and any centre that came out wrong has
  its two substituents swapped. A Fischer drawn from a wrong assignment is a
  convincing picture of the wrong sugar.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolTransforms import SetDihedralDeg

from .chemistry import ChemError, parse

FONT_STACK = "DejaVu Sans, Segoe UI, Arial, sans-serif"


@dataclass
class Projection:
    svg: str
    width: int
    height: int
    notes: list[str] = field(default_factory=list)


def _text(x: float, y: float, text: str, size: int, anchor: str = "middle") -> str:
    return (f'<text x="{x:.1f}" y="{y + size * 0.35:.1f}" font-family="{FONT_STACK}" '
            f'font-size="{size}" fill="#000000" text-anchor="{anchor}">{text}</text>')


def _line(x1: float, y1: float, x2: float, y2: float, width: float) -> str:
    return (f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
            f'stroke="#000000" stroke-width="{width:.2f}"/>')


def _svg(parts: list[str], width: float, height: float) -> str:
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0f}" '
            f'height="{height:.0f}" viewBox="0 0 {width:.0f} {height:.0f}">'
            + "".join(parts) + "</svg>")


# ------------------------------------------------------------------ Newman ---
def draw_newman(smiles: str, bond: list[int] | None = None, dihedral: float = 60.0,
                bond_px: int = 62, font_px: int = 18) -> Projection:
    mol = parse(smiles)
    if mol.GetNumBonds() == 0:
        raise ChemError("Newman projection needs at least one bond")
    work = Chem.AddHs(mol)
    front, back = _view_bond(mol, bond)

    status = AllChem.EmbedMolecule(work, randomSeed=0xF00D)
    if status != 0:
        AllChem.EmbedMolecule(work, useRandomCoords=True, randomSeed=0xF00D)

    conf = work.GetConformer()
    notes: list[str] = []
    _set_dihedral(work, front, back, dihedral, notes)

    axis = _sub(conf.GetAtomPosition(back), conf.GetAtomPosition(front))
    length = math.hypot(*axis) or 1.0
    axis = (axis[0] / length, axis[1] / length, axis[2] / length)
    e1 = _perpendicular(axis)
    e2 = _cross(axis, e1)

    front_dirs = _project(work, conf, front, back, axis, e1, e2)
    back_dirs = _project(work, conf, back, front, axis, e1, e2)

    radius = bond_px * 1.05
    outer = bond_px * 1.9
    size = int(outer * 2 + font_px * 6)
    cx = cy = size / 2

    parts: list[str] = [f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{radius:.1f}" '
                        f'fill="none" stroke="#000000" stroke-width="1.6"/>']

    for angle, atom in back_dirs:
        # In an eclipsed conformation a back bond sits exactly on a front one.
        # Nudge it sideways so both stay visible, keeping the true angle.
        nearest = min(abs(((angle - fa + math.pi) % (2 * math.pi)) - math.pi)
                      for fa, _ in front_dirs) if front_dirs else math.pi
        shift = 7.0 if nearest < math.radians(9) else 0.0
        ox, oy = -math.sin(angle) * shift, math.cos(angle) * shift
        x1 = cx + radius * math.cos(angle) + ox
        y1 = cy + radius * math.sin(angle) + oy
        x2 = cx + outer * math.cos(angle) + ox
        y2 = cy + outer * math.sin(angle) + oy
        parts.append(_line(x1, y1, x2, y2, 1.9))
        label = _atom_label(atom)
        if label:
            lx = cx + (outer + font_px * 0.8) * math.cos(angle) + ox
            ly = cy + (outer + font_px * 0.8) * math.sin(angle) + oy
            parts.append(_text(lx, ly, label, font_px, anchor="middle"))

    for angle, atom in front_dirs:
        x2, y2 = cx + radius * math.cos(angle), cy + radius * math.sin(angle)
        parts.append(_line(cx, cy, x2, y2, 2.3))
        label = _atom_label(atom)
        if label:
            # Front labels sit further out than the circle but closer than the
            # back labels, and shift the other way when eclipsed, so the two
            # never share a radius or a ray.
            nearest = min(abs(((angle - ba + math.pi) % (2 * math.pi)) - math.pi)
                          for ba, _ in back_dirs) if back_dirs else math.pi
            lateral = -7.0 if nearest < math.radians(9) else 0.0
            px, py = -math.sin(angle) * lateral, math.cos(angle) * lateral
            lx = cx + (radius + font_px * 1.7) * math.cos(angle) + px
            ly = cy + (radius + font_px * 1.7) * math.sin(angle) + py
            parts.append(_text(lx, ly, label, font_px, anchor="middle"))

    return Projection(svg=_svg(parts, size, size), width=size, height=size,
                      notes=notes)


def _atom_label(atom: Chem.Atom) -> str:
    """Label substituents in a Newman projection.

    A bare stick is standard skeletal shorthand, but in a Newman projection the
    whole point is reading what sits on each carbon, so a methyl is written out.
    """
    symbol = atom.GetSymbol()
    if symbol == "C":
        # Hydrogens are explicit after AddHs, so count neighbours, not
        # GetTotalNumHs, which is zero for every atom in that molecule.
        hydrogens = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() == 1)
        if hydrogens == 0:
            return ""
        return "CH" + (str(hydrogens) if hydrogens > 1 else "")
    if symbol in ("O", "S", "N"):
        return symbol
    return symbol


def _view_bond(mol: Chem.Mol, bond: list[int] | None) -> tuple[int, int]:
    if bond:
        if len(bond) != 2:
            raise ChemError("newman 'bond' needs two atom indices")
        i, j = int(bond[0]), int(bond[1])
        if not (0 <= i < mol.GetNumAtoms() and 0 <= j < mol.GetNumAtoms()):
            raise ChemError("newman 'bond' refers to an atom that does not exist")
        if mol.GetBondBetweenAtoms(i, j) is None:
            raise ChemError("newman 'bond' refers to atoms that are not bonded")
        return i, j

    # Degrees count hydrogens, so ethane qualifies: three substituents per carbon.
    counted = Chem.AddHs(mol)
    best = None
    for atom in counted.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        for neighbour in atom.GetNeighbors():
            if neighbour.GetAtomicNum() != 6:
                continue
            score = atom.GetDegree() + neighbour.GetDegree()
            if best is None or score > best[0]:
                best = (score, atom.GetIdx(), neighbour.GetIdx())
    if best is None:
        raise ChemError("Newman projection needs a bond carrying substituents")
    return best[1], best[2]


def _set_dihedral(mol, front, back, angle, notes) -> None:
    conf = mol.GetConformer()
    fa = [n.GetIdx() for n in mol.GetAtomWithIdx(front).GetNeighbors() if n.GetIdx() != back]
    ba = [n.GetIdx() for n in mol.GetAtomWithIdx(back).GetNeighbors() if n.GetIdx() != front]
    if not fa or not ba:
        notes.append("no substituent on one side, dihedral left as embedded")
        return
    try:
        SetDihedralDeg(conf, fa[0], front, back, ba[0], float(angle))
    except Exception:  # noqa: BLE001
        notes.append("could not set the requested dihedral")


def _project(mol, conf, centre, other, axis, e1, e2):
    origin = conf.GetAtomPosition(centre)
    out = []
    for neighbour in mol.GetAtomWithIdx(centre).GetNeighbors():
        if neighbour.GetIdx() == other:
            continue
        p = conf.GetAtomPosition(neighbour.GetIdx())
        v = (p.x - origin.x, p.y - origin.y, p.z - origin.z)
        dot = v[0] * axis[0] + v[1] * axis[1] + v[2] * axis[2]
        perp = (v[0] - dot * axis[0], v[1] - dot * axis[1], v[2] - dot * axis[2])
        a = perp[0] * e1[0] + perp[1] * e1[1] + perp[2] * e1[2]
        b = perp[0] * e2[0] + perp[1] * e2[1] + perp[2] * e2[2]
        out.append((math.atan2(b, a), neighbour))
    out.sort(key=lambda item: item[0])
    return out


def _sub(a, b):
    return (a.x - b.x, a.y - b.y, a.z - b.z)


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def _perpendicular(axis):
    helper = (0.0, 0.0, 1.0) if abs(axis[2]) < 0.9 else (1.0, 0.0, 0.0)
    e1 = _cross(axis, helper)
    n = math.hypot(*e1) or 1.0
    return (e1[0] / n, e1[1] / n, e1[2] / n)


# ----------------------------------------------------------------- Fischer ---
def draw_fischer(smiles: str, bond_px: int = 58, font_px: int = 17) -> Projection:
    mol = parse(smiles)
    if mol.GetRingInfo().NumRings() > 0:
        raise ChemError("Fischer projections are for open chains, not rings")
    work = Chem.AddHs(mol)
    chain = _longest_carbon_chain(work)
    if len(chain) < 3:
        raise ChemError("Fischer projection needs a chain of at least three carbons")
    chain = _orient_chain(work, chain)

    notes: list[str] = []
    coords, sides, ends = _fischer_layout(work, chain, bond_px, notes)
    _verify_and_fix(work, chain, coords, sides, bond_px, notes)
    return _draw_fischer(work, chain, coords, sides, ends, bond_px, font_px, notes)


def _orient_chain(mol: Chem.Mol, chain: list[int]) -> list[int]:
    """Fischer convention puts the most oxidised end at the top. Without this,
    glucose comes out upside down, which is the same molecule rotated 180° and
    still wrong for a textbook."""

    def oxidation(idx: int) -> float:
        score = 0.0
        for bond in mol.GetAtomWithIdx(idx).GetBonds():
            other = bond.GetOtherAtom(mol.GetAtomWithIdx(idx))
            if other.GetAtomicNum() == 8:
                score += bond.GetBondTypeAsDouble()
        return score

    if oxidation(chain[-1]) > oxidation(chain[0]):
        return list(reversed(chain))
    return chain


def _longest_carbon_chain(mol: Chem.Mol) -> list[int]:
    carbons = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6]
    best: list[int] = []

    def walk(path: list[int], seen: set[int]) -> None:
        nonlocal best
        if len(path) > len(best):
            best = list(path)
        for neighbour in mol.GetAtomWithIdx(path[-1]).GetNeighbors():
            idx = neighbour.GetIdx()
            if neighbour.GetAtomicNum() == 6 and idx not in seen:
                seen.add(idx)
                path.append(idx)
                walk(path, seen)
                path.pop()
                seen.remove(idx)

    for start in carbons:
        walk([start], {start})
    return best


def _fischer_layout(mol, chain, bond_px, notes):
    # Vertical bonds point away from the viewer and horizontal ones toward it.
    # The sign was calibrated against L-alanine and D-glucose: with horizontal
    # substituents at z = -1, RDKit reads the arrangement back as the stated
    # configuration, which is what makes this a Fischer projection rather than
    # a mirror image of one.
    VERTICAL_Z = 1.0
    HORIZONTAL_Z = -1.0

    coords: dict[int, tuple[float, float, float]] = {
        idx: (0.0, position * bond_px, VERTICAL_Z)
        for position, idx in enumerate(chain)
    }
    sides: dict[int, list[tuple[int, int]]] = {}
    ends: dict[int, str] = {}

    for position, idx in enumerate(chain):
        atom = mol.GetAtomWithIdx(idx)
        up = chain[position - 1] if position > 0 else None
        down = chain[position + 1] if position + 1 < len(chain) else None
        others = [n.GetIdx() for n in atom.GetNeighbors()
                  if n.GetIdx() not in (up, down)]
        chiral = atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED

        if position in (0, len(chain) - 1) and not chiral:
            ends[idx] = _group_label(mol, idx, up if up is not None else down)
            continue

        if len(others) == 2:
            a, b = others
            coords[a] = (-bond_px, position * bond_px, HORIZONTAL_Z)
            coords[b] = (bond_px, position * bond_px, HORIZONTAL_Z)
            sides[idx] = [(a, -1), (b, 1)]
        elif len(others) == 1:
            a = others[0]
            coords[a] = (-bond_px, position * bond_px, HORIZONTAL_Z)
            sides[idx] = [(a, -1)]
        else:
            notes.append(f"chain atom {idx} has no substituent to place")

    return coords, sides, ends


def _verify_and_fix(mol, chain, coords, sides, bond_px, notes) -> None:
    centres = [idx for idx in chain
               if mol.GetAtomWithIdx(idx).GetChiralTag()
               != Chem.ChiralType.CHI_UNSPECIFIED]
    if not centres:
        notes.append("no stereocentres, so the projection shows no configuration")
        return

    for _attempt in range(len(centres) + 2):
        conf = Chem.Conformer(mol.GetNumAtoms())
        conf.Set3D(True)
        for idx, (x, y, z) in coords.items():
            conf.SetAtomPosition(idx, (x, y, z * bond_px * 0.45))
        test = Chem.Mol(mol)
        test.RemoveAllConformers()
        test.AddConformer(conf, assignId=True)
        Chem.AssignAtomChiralTagsFromStructure(test, confId=0, replaceExistingTags=True)

        wrong = [idx for idx in centres
                 if test.GetAtomWithIdx(idx).GetChiralTag()
                 != mol.GetAtomWithIdx(idx).GetChiralTag()]
        if not wrong:
            return
        for idx in wrong:
            _swap_sides(coords, sides.get(idx, []))

    notes.append("could not resolve the Fischer arrangement for every centre")


def _swap_sides(coords, sides) -> None:
    """Mirror a centre's two horizontal substituents. Which atoms are the sides
    is recorded at layout time, not inferred from the coordinate sign, because
    the sign belongs to the projection convention."""
    for atom_idx, _sign in sides:
        x, y, z = coords[atom_idx]
        coords[atom_idx] = (-x, y, z)


def _group_label(mol, idx, parent) -> str:
    """Condensed label for a substituent, the way a Fischer projection writes
    it: CH3, OH, NH2, COOH, CH2OH."""
    atom = mol.GetAtomWithIdx(idx)
    symbol = atom.GetSymbol()
    hydrogens = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() == 1)
    heavy = [n for n in atom.GetNeighbors()
             if n.GetAtomicNum() > 1 and n.GetIdx() != parent]

    if symbol == "H":
        return "H"
    if symbol == "O":
        return "OH" if hydrogens else "O"
    if symbol == "N":
        if hydrogens > 1:
            return f"NH{hydrogens}"
        return "NH" if hydrogens == 1 else "N"
    if symbol != "C":
        return symbol + (f"H{hydrogens}" if hydrogens else "")

    if not heavy:
        return f"CH{hydrogens}" if hydrogens > 1 else ("CH" if hydrogens == 1 else "C")

    oxygens = [n for n in heavy if n.GetAtomicNum() == 8]
    if len(heavy) == 1 and len(oxygens) == 1:
        bond = mol.GetBondBetweenAtoms(idx, oxygens[0].GetIdx())
        double = bond is not None and bond.GetBondTypeAsDouble() == 2.0
        if double:
            return "CHO" if hydrogens == 1 else "C=O"
        return f"CH{hydrogens}OH" if hydrogens else "OH"
    if len(heavy) == 2 and len(oxygens) == 2:
        return "COOH"
    return symbol


def _draw_fischer(mol, chain, coords, sides, ends, bond_px, font_px, notes):
    xs = [c[0] for c in coords.values()]
    ys = [c[1] for c in coords.values()]
    margin = font_px * 6
    min_x, max_x = min(xs) - margin, max(xs) + margin
    min_y, max_y = min(ys) - font_px * 2, max(ys) + font_px * 2
    width = max(1.0, max_x - min_x)
    height = max(1.0, max_y - min_y)

    def X(x):
        return x - min_x

    def Y(y):
        return y - min_y

    parts: list[str] = []
    gap = font_px * 1.0
    for a, b in zip(chain, chain[1:]):
        if a not in coords or b not in coords:
            continue
        x1, y1 = coords[a][0], coords[a][1]
        x2, y2 = coords[b][0], coords[b][1]
        step = 1.0 if y2 >= y1 else -1.0
        if a in ends:
            y1 += step * gap
        if b in ends:
            y2 -= step * gap
        parts.append(_line(X(x1), Y(y1), X(x2), Y(y2), 2.0))

    for centre, items in sides.items():
        cx, cy = coords[centre][0], coords[centre][1]
        for atom_idx, sign in items:
            ax, ay, _ = coords[atom_idx]
            parts.append(_line(X(cx), Y(cy), X(ax), Y(ay), 2.0))
            label = _group_label(mol, atom_idx, centre)
            if ax < cx:
                parts.append(_text(X(ax) - 8, Y(ay), label, font_px, anchor="end"))
            else:
                parts.append(_text(X(ax) + 8, Y(ay), label, font_px, anchor="start"))

    for idx, label in ends.items():
        x, y, _ = coords[idx]
        parts.append(_text(X(x), Y(y), label, font_px))

    return Projection(svg=_svg(parts, width, height), width=int(width),
                      height=int(height), notes=notes)








