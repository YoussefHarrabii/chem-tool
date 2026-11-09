"""Newman and Fischer projections.

The point of these tests is stereochemistry: a projection that is mirrored, or
one where a centre was placed wrong, is a convincing picture of the wrong
molecule. The layout is checked by reading it back, not by comparing pixels.
"""

from __future__ import annotations

import pytest
from rdkit import Chem

from chemdraft.chemistry import ChemError
from chemdraft.projection import draw_fischer, draw_newman
from chemdraft.scene import SpecError, render_scene


def scene(items):
    return {"style": {"reference": False}, "items": items}


# ------------------------------------------------------------------ Newman ---
def test_newman_renders():
    projection = draw_newman("CCCC", bond=[1, 2], dihedral=60)
    assert projection.width > 0 and projection.height > 0
    assert "<svg" in projection.svg


def test_newman_picks_a_bond_when_none_given():
    projection = draw_newman("CCCC")
    assert projection.width > 0


def test_newman_ethane():
    """Ethane has no bond with degree 2 on both ends until hydrogens count."""
    projection = draw_newman("CC")
    assert projection.width > 0


def test_newman_conformations_differ():
    staggered = draw_newman("CCCC", bond=[1, 2], dihedral=60).svg
    eclipsed = draw_newman("CCCC", bond=[1, 2], dihedral=0).svg
    assert staggered != eclipsed


def test_newman_rejects_a_bond_that_is_not_there():
    with pytest.raises(ChemError):
        draw_newman("CCCC", bond=[0, 3])


def test_newman_rejects_an_out_of_range_atom():
    with pytest.raises(ChemError):
        draw_newman("CCCC", bond=[0, 99])


def _methyl_separation(smiles: str, bond, dihedral: float) -> float:
    """Angle between the two methyl groups as drawn, in degrees."""
    import math

    from rdkit.Chem import AllChem
    from rdkit.Chem.rdMolTransforms import SetDihedralDeg

    from chemdraft.projection import (_cross, _perpendicular, _project, _sub,
                                      _view_bond)

    mol = Chem.MolFromSmiles(smiles)
    work = Chem.AddHs(mol)
    front, back = _view_bond(mol, bond)
    AllChem.EmbedMolecule(work, randomSeed=0xF00D)
    conf = work.GetConformer()
    fa = [n.GetIdx() for n in work.GetAtomWithIdx(front).GetNeighbors()
          if n.GetIdx() != back]
    ba = [n.GetIdx() for n in work.GetAtomWithIdx(back).GetNeighbors()
          if n.GetIdx() != front]
    SetDihedralDeg(conf, fa[0], front, back, ba[0], dihedral)

    axis = _sub(conf.GetAtomPosition(back), conf.GetAtomPosition(front))
    length = math.hypot(*axis)
    axis = (axis[0] / length, axis[1] / length, axis[2] / length)
    e1 = _perpendicular(axis)
    e2 = _cross(axis, e1)

    def methyl_angle(centre, other):
        for angle, atom in _project(work, conf, centre, other, axis, e1, e2):
            if atom.GetAtomicNum() == 6:
                return math.degrees(angle) % 360
        raise AssertionError("no methyl found")

    a = methyl_angle(front, back)
    b = methyl_angle(back, front)
    return abs((a - b + 180) % 360 - 180)


@pytest.mark.parametrize("requested", [180, 120, 60, 0])
def test_newman_dihedral_is_what_the_drawing_shows(requested):
    """The whole point of the projection: the angle a reader sees must be the
    angle that was asked for. Gemini read panel 1 (anti) and panel 3 (eclipsed)
    as both 'staggered 60', so this is measured, not eyeballed."""
    assert abs(_methyl_separation("CCCC", [1, 2], requested) - requested) < 1.0


# ----------------------------------------------------------------- Fischer ---
def test_fischer_renders_alanine():
    projection = draw_fischer("C[C@H](N)C(=O)O")
    assert projection.width > 0
    assert not projection.notes


def test_fischer_rejects_a_ring():
    with pytest.raises(ChemError):
        draw_fischer("C1CCCCC1")


def test_fischer_reports_no_stereocentres():
    projection = draw_fischer("CCCC")
    assert any("no stereocentres" in note for note in projection.notes)


def test_fischer_distinguishes_enantiomers():
    """L and D alanine must not produce the same drawing, and neither may fail
    the arrangement check."""
    left = draw_fischer("C[C@H](N)C(=O)O")
    right = draw_fischer("C[C@@H](N)C(=O)O")
    assert left.svg != right.svg
    assert not left.notes and not right.notes


def _read_back(smiles: str) -> list[str]:
    """Build the Fischer arrangement and ask RDKit what configuration it sees."""
    from chemdraft.projection import (_fischer_layout, _longest_carbon_chain,
                                      _orient_chain, _verify_and_fix)

    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    chain = _orient_chain(mol, _longest_carbon_chain(mol))
    notes: list[str] = []
    coords, sides, _ends = _fischer_layout(mol, chain, 50, notes)
    _verify_and_fix(mol, chain, coords, sides, 50, notes)

    conf = Chem.Conformer(mol.GetNumAtoms())
    conf.Set3D(True)
    for idx, (x, y, z) in coords.items():
        conf.SetAtomPosition(idx, (x, y, z * 50 * 0.45))
    test = Chem.Mol(mol)
    test.RemoveAllConformers()
    test.AddConformer(conf, assignId=True)
    Chem.AssignAtomChiralTagsFromStructure(test, confId=0, replaceExistingTags=True)
    Chem.AssignStereochemistry(test, cleanIt=True, force=True)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    found = Chem.FindMolChiralCenters(test, includeUnassigned=False,
                                      useLegacyImplementation=False)
    return sorted(label for _, label in found)


def test_fischer_alanine_reads_back_as_s():
    assert _read_back("C[C@H](N)C(=O)O") == ["S"]


def test_fischer_glucose_reads_back_as_d_glucose():
    """D-glucose open chain is 2R,3S,4R,5R."""
    labels = _read_back("O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO")
    assert labels == ["R", "R", "R", "S"]


# ------------------------------------------------------------------- scene ---
def test_newman_and_fischer_in_a_scene():
    result = render_scene(scene([
        {"type": "row", "items": [
            {"type": "newman", "smiles": "CCCC", "bond": [1, 2], "dihedral": 60},
            {"type": "fischer", "smiles": "C[C@H](N)C(=O)O"},
        ]}
    ]))
    assert result.png[:4] == b"\x89PNG"


def test_projection_needs_smiles():
    with pytest.raises(SpecError):
        render_scene(scene([{"type": "fischer"}]))
