"""Choose a layout instead of accepting the first one.

2D layout is a search problem with no single answer. CoordGen is usually best,
but not always: on bridged bicyclics it can place a substituent almost on top
of its parent atom, which draws a stub instead of a bond. RDKit's older
depiction with ring templates handles those better and worse elsewhere.

So: generate candidates, score them on properties that make a drawing readable,
and keep the best. The scoring contains no knowledge of any molecule, and the
same four candidates are tried for everything.
"""

from __future__ import annotations

import itertools
import math
import statistics

from rdkit import Chem
from rdkit.Chem import rdDepictor

# Candidates: (prefer CoordGen, use ring templates)
CANDIDATES = [
    (True, False),
    (True, True),
    (False, False),
    (False, True),
]

# Layout is skipped for very large molecules: quadratic scoring costs more than
# it buys, and CoordGen is the sensible default there.
SCORING_LIMIT = 220


def _points(mol: Chem.Mol) -> list[tuple[float, float]]:
    conf = mol.GetConformer()
    return [(conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y)
            for i in range(mol.GetNumAtoms())]


def _orient(a, b, c) -> float:
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def _crosses(p1, p2, p3, p4) -> bool:
    d1, d2 = _orient(p3, p4, p1), _orient(p3, p4, p2)
    d3, d4 = _orient(p1, p2, p3), _orient(p1, p2, p4)
    return ((d1 > 0) != (d2 > 0)) and ((d3 > 0) != (d4 > 0))


def _point_segment_distance(p, a, b) -> float:
    ax, ay = a
    bx, by = b
    dx, dy = bx - ax, by - ay
    length2 = dx * dx + dy * dy
    if length2 == 0:
        return math.dist(p, a)
    t = max(0.0, min(1.0, ((p[0] - ax) * dx + (p[1] - ay) * dy) / length2))
    return math.dist(p, (ax + t * dx, ay + t * dy))


def score_layout(mol: Chem.Mol) -> tuple[int, float, int, float]:
    """Lower is better. Ordered by how badly each problem ruins a figure:
    bonds crossing, atoms colliding, atoms sitting on a bond that is not
    theirs, then uneven bond lengths."""
    pts = _points(mol)
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]

    crossings = 0
    for (a1, b1), (a2, b2) in itertools.combinations(bonds, 2):
        if {a1, b1} & {a2, b2}:
            continue
        if _crosses(pts[a1], pts[b1], pts[a2], pts[b2]):
            crossings += 1

    lengths = [math.dist(pts[a], pts[b]) for a, b in bonds] or [1.0]
    median = statistics.median(lengths)

    collisions = 0.0
    bonded = {frozenset(b) for b in bonds}
    for i, j in itertools.combinations(range(len(pts)), 2):
        if frozenset((i, j)) in bonded:
            continue
        distance = math.dist(pts[i], pts[j])
        if distance < 0.6 * median:
            # Weight by how badly they overlap, not just by how many pairs do.
            # Two atoms a hair apart must never beat one pair slightly tight.
            collisions += (0.6 * median - distance) / median

    on_bond = 0
    for atom in range(len(pts)):
        for a, b in bonds:
            if atom in (a, b):
                continue
            if _point_segment_distance(pts[atom], pts[a], pts[b]) < 0.18 * median:
                on_bond += 1

    spread = statistics.pstdev(lengths) / median if median else 0.0
    return (crossings, collisions, on_bond, round(spread, 3))


def compute_layout(mol: Chem.Mol) -> tuple[Chem.Mol, str, tuple]:
    """Return the best of the candidate layouts, its name, and its score."""
    if mol.GetNumBonds() > SCORING_LIMIT:
        work = Chem.Mol(mol)
        rdDepictor.SetPreferCoordGen(True)
        rdDepictor.Compute2DCoords(work)
        return work, "coordgen (large molecule, not scored)", (0, 0, 0, 0.0)

    best = None
    for coordgen, ring_templates in CANDIDATES:
        work = Chem.Mol(mol)
        work.RemoveAllConformers()
        rdDepictor.SetPreferCoordGen(coordgen)
        try:
            rdDepictor.Compute2DCoords(work, useRingTemplates=ring_templates)
        except Exception:  # noqa: BLE001 - a candidate that fails is just not offered
            continue
        name = ("coordgen" if coordgen else "rdkit")
        name += " + ring templates" if ring_templates else ""
        result = score_layout(work)
        if best is None or result < best[2]:
            best = (work, name, result)

    rdDepictor.SetPreferCoordGen(True)  # leave global state as we found it

    if best is None:
        work = Chem.Mol(mol)
        rdDepictor.Compute2DCoords(work)
        return work, "coordgen (all candidates failed)", (0, 0, 0, 0.0)
    return best
