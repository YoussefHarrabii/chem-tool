"""The faithfulness gate.

Every figure is checked against the molecule it was built from before a single
byte is written. If the check cannot prove the drawing shows the parsed
molecule, the tool refuses to produce a file. A missing figure is safe; a wrong
figure teaches the wrong thing.
"""

from __future__ import annotations

import itertools

import re
from dataclasses import dataclass, field

import numpy as np
from rdkit import Chem

from .chemistry import Description, describe, needs_label

ATOM_CLASS = re.compile(r"class='atom-(\d+)'")
BOND_CLASS = re.compile(r"class='bond-(\d+)")


@dataclass
class Report:
    problems: list[str]
    checks_run: list[str]
    notes: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.problems


def _svg_index(svg: str) -> tuple[set[int], set[int]]:
    """Atom labels and bonds that actually appear in the drawing."""
    labelled, bonds = set(), set()
    for m in ATOM_CLASS.finditer(svg):
        labelled.add(int(m.group(1)))
    for m in BOND_CLASS.finditer(svg):
        bonds.add(int(m.group(1)))
    return labelled, bonds


def _side(point, origin, axis) -> float:
    v = (point[0] - origin[0], point[1] - origin[1])
    return axis[0] * v[1] - axis[1] * v[0]


def verify(
    parsed: Chem.Mol,
    prepared: Chem.Mol,
    drawer,
    svg: str,
    width: int,
    height: int,
    label_all: bool,
    expected: Description,
) -> Report:
    problems: list[str] = []
    checks: list[str] = []
    notes: list[str] = []

    # --- 1. the drawing's molecule is the parsed molecule ----------------
    got = describe(prepared)
    for field_name in ("formula", "charge", "radicals"):
        if getattr(got, field_name) != getattr(expected, field_name):
            problems.append(
                f"{field_name} changed during preparation: "
                f"input {getattr(expected, field_name)!r}, "
                f"drawn {getattr(got, field_name)!r}"
            )
    checks.append("formula, charge and radical count survive preparation")

    if len(got.stereocentres) != len(expected.stereocentres):
        problems.append(
            f"stereocentre count changed: input {expected.stereocentres}, "
            f"drawn {got.stereocentres}"
        )
    if got.double_bond_stereo != expected.double_bond_stereo:
        problems.append(
            f"double bond stereo count changed: input "
            f"{expected.double_bond_stereo}, drawn {got.double_bond_stereo}"
        )
    checks.append("stereochemistry survives preparation")

    # --- 2. every bond is actually drawn ---------------------------------
    labelled, bonds_drawn = _svg_index(svg)
    missing_bonds = [b.GetIdx() for b in prepared.GetBonds()
                     if b.GetIdx() not in bonds_drawn]
    if missing_bonds:
        problems.append(f"bonds missing from the drawing: {missing_bonds}")
    checks.append("every bond appears in the drawing")

    # --- 3. every atom that must be named is named -----------------------
    missing_labels = []
    for atom in prepared.GetAtoms():
        idx = atom.GetIdx()
        must = label_all or needs_label(atom)
        if must and idx not in labelled:
            missing_labels.append(
                f"atom {idx} ({atom.GetSymbol()}, charge {atom.GetFormalCharge()}, "
                f"H{atom.GetTotalNumHs()})"
            )
    if missing_labels:
        problems.append(
            "atoms that must show a symbol but do not: "
            + ", ".join(missing_labels)
        )
    checks.append("every atom needing a symbol is labelled in the drawing")

    # --- 4. assigned stereocentres have a wedge or hash -------------------
    wedged = {
        b.GetBeginAtomIdx() for b in prepared.GetBonds()
        if b.GetBondDir() in (Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH)
    } | {
        b.GetEndAtomIdx() for b in prepared.GetBonds()
        if b.GetBondDir() in (Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH)
    }
    for idx, _cfg in Chem.FindMolChiralCenters(
        prepared, includeUnassigned=False, useLegacyImplementation=False
    ):
        if idx not in wedged:
            problems.append(
                f"stereocentre at atom {idx} has no wedge or hash, so the "
                f"drawing implies the wrong configuration"
            )
    checks.append("assigned stereocentres are drawn with a wedge or hash")

    # --- 4c. the wedges encode the stated configuration -------------------
    if expected.stereocentres:
        try:
            replayed = Chem.Mol(prepared)
            for atom in replayed.GetAtoms():
                atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
                if atom.HasProp("_CIPCode"):
                    atom.ClearProp("_CIPCode")
            Chem.AssignChiralTypesFromBondDirs(replayed)
            Chem.AssignStereochemistry(replayed, cleanIt=True, force=True)
            Chem.AssignStereochemistry(prepared, cleanIt=True, force=True)
            before = sorted(label for _, label in Chem.FindMolChiralCenters(
                prepared, includeUnassigned=False, useLegacyImplementation=False))
            after = sorted(label for _, label in Chem.FindMolChiralCenters(
                replayed, includeUnassigned=False, useLegacyImplementation=False))
            if before != after:
                problems.append(
                    f"the drawn wedges encode {after}, but the molecule is "
                    f"{before}, so the figure shows the wrong stereoisomer"
                )
            checks.append("wedges reproduce the stated configuration")
        except Exception as exc:  # noqa: BLE001 - older RDKit may lack the call
            notes.append(f"could not read the wedges back: {exc}")

    # --- 4b. say when the layout chose a geometry the input left open -----
    try:
        for info in Chem.FindPotentialStereo(prepared):
            if info.specified != Chem.StereoSpecified.Unspecified:
                continue
            if info.type == Chem.StereoType.Bond_Double:
                notes.append(
                    f"double bond {info.centeredOn} had no specified geometry, "
                    f"so the drawing shows one particular isomer"
                )
            elif info.type == Chem.StereoType.Atom_Tetrahedral:
                notes.append(
                    f"stereocentre at atom {info.centeredOn} is unspecified and "
                    f"is drawn flat"
                )
    except Exception:  # noqa: BLE001 - older RDKit lacks FindPotentialStereo
        notes.append("could not check for unspecified stereochemistry")

    # --- 5. E/Z drawn on the correct side --------------------------------
    conf = prepared.GetConformer()
    for bond in prepared.GetBonds():
        stereo = bond.GetStereo()
        if stereo == Chem.BondStereo.STEREONONE:
            continue
        if stereo == Chem.BondStereo.STEREOANY:
            problems.append(
                f"double bond {bond.GetIdx()} has unspecified stereo but is "
                f"drawn as if it were definite"
            )
            continue
        want_cis = stereo in (Chem.BondStereo.STEREOCIS, Chem.BondStereo.STEREOZ)
        want_trans = stereo in (Chem.BondStereo.STEREOTRANS, Chem.BondStereo.STEREOE)
        if not (want_cis or want_trans):
            continue
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        # A cis/trans descriptor refers to specific substituents. RDKit records
        # which ones; guessing a neighbour reports false alarms on any alkene
        # carrying two groups on one end.
        refs = list(bond.GetStereoAtoms())
        if len(refs) != 2:
            ni = [n.GetIdx() for n in prepared.GetAtomWithIdx(i).GetNeighbors()
                  if n.GetIdx() != j]
            nj = [n.GetIdx() for n in prepared.GetAtomWithIdx(j).GetNeighbors()
                  if n.GetIdx() != i]
            if len(ni) == 1 and len(nj) == 1:
                refs = [ni[0], nj[0]]
            else:
                notes.append(
                    f"double bond {bond.GetIdx()}: stereo not checked, no "
                    f"reference substituents recorded"
                )
                continue
        pi, pj = conf.GetAtomPosition(i), conf.GetAtomPosition(j)
        axis = (pj.x - pi.x, pj.y - pi.y)
        if abs(axis[0]) < 1e-6 and abs(axis[1]) < 1e-6:
            problems.append(f"double bond {bond.GetIdx()} has degenerate geometry")
            continue
        a = conf.GetAtomPosition(refs[0])
        b = conf.GetAtomPosition(refs[1])
        sign = _side((a.x, a.y), (pi.x, pi.y), axis) * _side(
            (b.x, b.y), (pi.x, pi.y), axis
        )
        drawn_cis = sign > 0
        if drawn_cis != want_cis:
            problems.append(
                f"double bond {bond.GetIdx()} is {'cis' if want_cis else 'trans'} "
                f"in the input but drawn {'cis' if drawn_cis else 'trans'}"
            )
    checks.append("double bond geometry matches the declared E/Z")

    # --- 6. nothing is drawn off the canvas ------------------------------
    outside = []
    for atom in prepared.GetAtoms():
        p = drawer.GetDrawCoords(atom.GetIdx())
        if not (0 <= p.x <= width and 0 <= p.y <= height):
            outside.append(f"atom {atom.GetIdx()} at ({p.x:.0f},{p.y:.0f})")
    if outside:
        problems.append("atoms drawn outside the canvas: " + ", ".join(outside))
    checks.append("every atom lies inside the canvas")

    # --- 7. no two atoms are drawn on top of each other ------------------
    pts = [(conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y)
           for i in range(prepared.GetNumAtoms())]
    lengths = []
    for bond in prepared.GetBonds():
        a = conf.GetAtomPosition(bond.GetBeginAtomIdx())
        b = conf.GetAtomPosition(bond.GetEndAtomIdx())
        lengths.append(((a.x - b.x) ** 2 + (a.y - b.y) ** 2) ** 0.5)
    median = sorted(lengths)[len(lengths) // 2] if lengths else 1.0
    bonded = {frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
              for b in prepared.GetBonds()}
    for i, j in itertools.combinations(range(len(pts)), 2):
        if frozenset((i, j)) in bonded:
            continue
        distance = ((pts[i][0] - pts[j][0]) ** 2
                    + (pts[i][1] - pts[j][1]) ** 2) ** 0.5
        if distance < 0.2 * median:
            problems.append(
                f"atoms {i} and {j} are drawn on top of each other, so the "
                f"figure cannot be read unambiguously"
            )
        elif distance < 0.45 * median:
            notes.append(
                f"atoms {i} and {j} are drawn close together "
                f"({distance / median:.2f} bond lengths)"
            )
    checks.append("no two atoms are drawn on top of each other")

    return Report(problems=problems, checks_run=checks, notes=notes)


