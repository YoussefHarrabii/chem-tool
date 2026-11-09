"""Tests for the faithfulness gate.

The point of these is not coverage. It is to prove the gate has teeth: each
test breaks a drawing on purpose and asserts the gate refuses it. A check that
has never been seen to fail is not a check.
"""

from __future__ import annotations

import numpy as np
import pytest
from rdkit import Chem

from chemdraft.chemistry import ChemError, describe, parse
from chemdraft.render import RenderRequest, render
from chemdraft.verify import verify


# ---------------------------------------------------------------- helpers ---
def drawn(smiles: str, mode: str = "skeletal"):
    return render(RenderRequest(smiles=smiles, mode=mode, target_width=600))


def recheck(result, **overrides):
    """Re-run the gate on an already drawn figure, optionally lying about the
    canvas or the SVG so the gate has something wrong to find."""
    args = dict(
        parsed=result.prepared,
        prepared=result.prepared,
        drawer=result.drawer,
        svg=result.svg,
        width=result.width,
        height=result.height,
        label_all=False,
        expected=result.description,
    )
    args.update(overrides)
    return verify(**args)


# ------------------------------------------------------ rejecting bad input ---
@pytest.mark.parametrize(
    "smiles",
    [
        "C(C)(C)(C)(C)C",   # pentavalent carbon
        "C1CC",             # unclosed ring
        "c1cccc1",          # invalid aromatic ring
        "CC(=O",            # unbalanced bracket
        "",                 # nothing at all
    ],
)
def test_bad_input_raises(smiles):
    with pytest.raises(ChemError):
        parse(smiles)


def test_no_file_written_on_failure(tmp_path, capsys):
    from chemdraft.cli import main

    out = tmp_path / "nope.png"
    code = main(["smiles", "C(C)(C)(C)(C)C", "-o", str(out)])
    assert code == 2
    assert not out.exists()


# ------------------------------------------------------- the gate has teeth ---
def test_gate_catches_flipped_double_bond():
    """Draw a trans alkene, then move a substituent to the other side. The gate
    must notice, because otherwise a trans compound could be published as cis."""
    r = drawn("F/C=C/F")
    assert r.report.ok

    conf = r.prepared.GetConformer()
    double = [b for b in r.prepared.GetBonds()
              if b.GetBondType() == Chem.BondType.DOUBLE][0]
    i, j = double.GetBeginAtomIdx(), double.GetEndAtomIdx()
    other = [n.GetIdx() for n in r.prepared.GetAtomWithIdx(j).GetNeighbors()
             if n.GetIdx() != i][0]

    # Reflect the substituent across the double bond axis. That is exactly the
    # cis/trans flip, and mirroring on an arbitrary axis would have done nothing.
    pi_at, pj_at, p_at = (conf.GetAtomPosition(k) for k in (i, j, other))
    pi = np.array([pi_at.x, pi_at.y])
    pj = np.array([pj_at.x, pj_at.y])
    p = np.array([p_at.x, p_at.y])
    u = pj - pi
    u = u / np.linalg.norm(u)
    v = p - pi
    perp = v - np.dot(v, u) * u
    assert np.linalg.norm(perp) > 1e-3, "substituent sits on the axis, no flip possible"
    flipped = pi + v - 2 * perp
    conf.SetAtomPosition(
        other, Chem.rdGeometry.Point3D(float(flipped[0]), float(flipped[1]), 0.0)
    )

    report = recheck(r)
    assert not report.ok
    assert any("drawn cis" in msg or "drawn trans" in msg for msg in report.problems)


def test_gate_catches_missing_bond_class():
    """If the SVG lost a bond, the drawing is not the molecule."""
    r = drawn("CCO")
    stripped = r.svg
    for n in range(r.prepared.GetNumBonds()):
        stripped = stripped.replace(f"bond-{n}", f"bondx-{n}")
    report = recheck(r, svg=stripped)
    assert not report.ok
    assert any("bonds missing" in msg for msg in report.problems)


def test_gate_catches_missing_label():
    """A nitrogen drawn as a bare vertex would be read as a carbon."""
    r = drawn("CCO")
    stripped = r.svg.replace("class='atom-2'", "class='nope-2'")
    report = recheck(r, svg=stripped)
    assert not report.ok
    assert any("must show a symbol" in msg for msg in report.problems)


def test_gate_catches_clipped_canvas():
    r = drawn("CCO")
    report = recheck(r, width=10, height=10)
    assert not report.ok
    assert any("outside the canvas" in msg for msg in report.problems)


def test_gate_catches_lost_stereocentre():
    """If preparation dropped the wedge, the drawing shows a flat carbon where
    the input had a definite configuration."""
    r = drawn("C[C@H](N)C(=O)O")
    assert r.report.ok
    for bond in r.prepared.GetBonds():
        bond.SetBondDir(Chem.BondDir.NONE)
    report = recheck(r)
    assert not report.ok
    assert any("no wedge or hash" in msg for msg in report.problems)


def test_gate_catches_wrong_wedge():
    """A wedge on the wrong face draws the other enantiomer. The gate must
    notice, because both drawings look equally plausible."""
    r = drawn("C[C@H](N)C(=O)O")
    assert r.report.ok

    flipped = False
    for bond in r.prepared.GetBonds():
        if bond.GetBondDir() == Chem.BondDir.BEGINWEDGE:
            bond.SetBondDir(Chem.BondDir.BEGINDASH)
            flipped = True
        elif bond.GetBondDir() == Chem.BondDir.BEGINDASH:
            bond.SetBondDir(Chem.BondDir.BEGINWEDGE)
            flipped = True
    assert flipped, "no wedge found to flip"

    report = recheck(r)
    assert not report.ok
    assert any("wrong stereoisomer" in problem for problem in report.problems)


# ------------------------------------------------------------- happy paths ---
@pytest.mark.parametrize(
    "smiles,formula",
    [
        ("CC(=O)O", "C2H4O2"),
        ("CC(=O)N", "C2H5NO"),
        ("c1ccccc1", "C6H6"),
        ("C[C@H](N)C(=O)O", "C3H7NO2"),
        ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "C8H10N4O2"),   # caffeine
        ("CC(=O)Oc1ccccc1C(=O)O", "C9H8O4"),            # aspirin
        ("[NH4+]", "H4N+"),
        ("CC(=O)[O-]", "C2H3O2-"),
    ],
)
def test_formula_is_derived_not_guessed(smiles, formula):
    assert describe(parse(smiles)).formula == formula


@pytest.mark.parametrize("smiles", ["F/C=C/F", "F/C=C\\F"])
def test_double_bond_stereo_is_kept(smiles):
    mol = parse(smiles)
    assert describe(mol).double_bond_stereo == 1, "stereo was dropped at parse time"
    assert parse(smiles) is not None


def test_stereocentre_survives_to_the_drawing():
    r = drawn("C[C@H](N)C(=O)O")
    assert r.description.stereocentres == ["1:S"]
    assert r.report.ok


def test_lewis_mode_labels_every_atom():
    r = drawn("C", mode="lewis")
    assert r.report.ok
    for atom in r.prepared.GetAtoms():
        assert atom.HasProp("atomLabel") or atom.GetSymbol() != "C"


def test_canonical_round_trip():
    for smiles in ["C[C@H](N)C(=O)O", "F/C=C/F", "c1ccc2ccccc2c1", "[O-][N+](=O)C"]:
        mol = parse(smiles)
        canonical = Chem.MolToSmiles(mol)
        assert Chem.MolToSmiles(Chem.MolFromSmiles(canonical)) == canonical
