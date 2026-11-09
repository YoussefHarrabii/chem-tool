"""Lone pairs, drawn as dots.

RDKit never draws lone pairs, so a carbonyl oxygen comes out with no dots and
looks wrong to anyone who has been taught Lewis structures.

The count is arithmetic, not a lookup table:

    non-bonding electrons = outer electrons of the element
                          - formal charge
                          - (bond orders + hydrogens)

The placement is arithmetic too: dots go in the widest free angular gap around
the atom, outside its label.
"""

from __future__ import annotations

import math

from rdkit import Chem

_PT = Chem.GetPeriodicTable()


def nonbonding_electrons(mol: Chem.Mol, idx: int) -> tuple[int, int]:
    """(pairs, unpaired) for one atom. Works on a kekulised copy so aromatic
    bonds are whole numbers."""
    atom = mol.GetAtomWithIdx(idx)
    order = 0.0
    for bond in atom.GetBonds():
        kind = bond.GetBondTypeAsDouble()
        order += 1.5 if kind == 1.5 else kind
    order += atom.GetTotalNumHs()
    outer = _PT.GetNOuterElecs(atom.GetAtomicNum())
    left = outer - atom.GetFormalCharge() - order
    if left <= 0:
        return 0, 0
    n = int(round(left))
    return n // 2, n % 2


def _angle_gaps(mol: Chem.Mol, idx: int, drawer) -> list[tuple[float, float]]:
    """Free angular gaps around an atom in screen space, widest first."""
    centre = drawer.GetDrawCoords(idx)
    angles = []
    for bond in mol.GetAtomWithIdx(idx).GetBonds():
        other = drawer.GetDrawCoords(bond.GetOtherAtomIdx(idx))
        angles.append(math.atan2(other.y - centre.y, other.x - centre.x))

    if not angles:
        return [(2 * math.pi, 0.0)]

    angles.sort()
    gaps = []
    for i, start in enumerate(angles):
        end = angles[(i + 1) % len(angles)]
        span = (end - start) % (2 * math.pi)
        if span == 0:
            span = 2 * math.pi
        gaps.append((span, start + span / 2))
    gaps.sort(key=lambda g: g[0], reverse=True)
    return gaps


def _slot_angles(gaps: list[tuple[float, float]], count: int) -> list[float]:
    """Pick `count` angles, spreading through the free gaps. An atom with a
    single neighbour has one gap wrapping the whole circle, so several pairs
    must be fanned out inside it or they land on top of each other."""
    if count <= 0:
        return []
    if count <= len(gaps):
        return [gaps[i][1] for i in range(count)]

    angles: list[float] = []
    remaining = count
    for index, (span, mid) in enumerate(gaps):
        if remaining == 0:
            break
        share = remaining if index == len(gaps) - 1 else 1
        share = min(share, remaining)
        step = min(0.85, span / (share + 1))
        for j in range(share):
            angles.append(mid + (j - (share - 1) / 2) * step)
        remaining -= share
    while len(angles) < count:
        angles.append(gaps[0][1])
    return angles


def dots_for(mol: Chem.Mol, drawer, font_px: float) -> list[str]:
    """SVG circles for every lone pair in the prepared molecule."""
    circles: list[str] = []
    radius = max(1.1, font_px * 0.055)
    offset = font_px * 0.72          # distance from the atom centre
    separation = font_px * 0.30      # distance between the two dots of a pair

    for atom in mol.GetAtoms():
        pairs, unpaired = nonbonding_electrons(mol, atom.GetIdx())
        if pairs == 0 and unpaired == 0:
            continue

        slots = [2] * pairs + ([1] if unpaired else [])
        gaps = _angle_gaps(mol, atom.GetIdx(), drawer)
        angles = _slot_angles(gaps, len(slots))
        centre = drawer.GetDrawCoords(atom.GetIdx())

        for dots, mid in zip(slots, angles):
            for d in range(dots):
                lateral = (d - (dots - 1) / 2) * separation
                x = centre.x + offset * math.cos(mid) - lateral * math.sin(mid)
                y = centre.y + offset * math.sin(mid) + lateral * math.cos(mid)
                circles.append(
                    f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{radius:.2f}" '
                    f'fill="#000000"/>'
                )
    return circles


def inject(svg: str, circles: list[str]) -> str:
    if not circles:
        return svg
    marker = "</svg>"
    at = svg.rfind(marker)
    if at < 0:
        return svg
    return svg[:at] + "\n".join(circles) + "\n" + svg[at:]
