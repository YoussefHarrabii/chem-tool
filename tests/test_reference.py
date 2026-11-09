"""Reference depictions.

The point of these tests is that the reference path is either exact or absent.
A wrong record must never be drawn, and a missing record must fall back to
automatic layout with the caller told about it.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import inchi

from chemdraft import reference
from chemdraft.render import RenderRequest, render

CAMPHOR = "CC12CCC(CC1=O)C2(C)C"
FIXTURE = Path(__file__).parent / "data" / "camphor-2d.sdf"


def key_of(smiles: str) -> str:
    return inchi.MolToInchiKey(Chem.MolFromSmiles(smiles))


@pytest.fixture
def cache(tmp_path, monkeypatch):
    monkeypatch.setattr(reference, "CACHE", tmp_path)
    return tmp_path


def install(cache: Path, key: str, sdf: str) -> None:
    (cache / f"{key}.sdf").write_text(sdf, encoding="utf-8")


def test_reference_is_used_when_it_matches(cache):
    install(cache, key_of(CAMPHOR), FIXTURE.read_text(encoding="utf-8"))
    result = render(RenderRequest(smiles=CAMPHOR, reference=True))
    assert result.report.ok
    assert result.depiction.startswith("reference depiction")

    # The drawn coordinates are the reference coordinates.
    ref = Chem.MolFromMolBlock(FIXTURE.read_text(encoding="utf-8"), removeHs=True)
    match = result.prepared.GetSubstructMatch(ref)
    assert len(match) == ref.GetNumAtoms()
    ref_conf, conf = ref.GetConformer(), result.prepared.GetConformer()
    for ref_idx, mol_idx in enumerate(match):
        a, b = ref_conf.GetAtomPosition(ref_idx), conf.GetAtomPosition(mol_idx)
        assert abs(a.x - b.x) < 1e-6 and abs(a.y - b.y) < 1e-6


def test_wrong_record_is_refused(cache):
    """A record for a different molecule must not be drawn."""
    wrong = Chem.MolToMolBlock(Chem.MolFromSmiles("CC(=O)O"))
    install(cache, key_of(CAMPHOR), wrong)
    result = render(RenderRequest(smiles=CAMPHOR, reference=True))
    assert result.depiction.startswith("automatic")
    assert any("does not match" in note for note in result.notes)


def test_missing_record_falls_back(cache, monkeypatch):
    monkeypatch.setattr(reference, "_download", lambda key: (None, "no reference available (offline)"))
    result = render(RenderRequest(smiles="CCO", reference=True))
    assert result.report.ok
    assert result.depiction.startswith("automatic")
    assert any("offline" in note for note in result.notes)


def test_polycyclic_skeleton_warns_when_automatic(cache, monkeypatch):
    monkeypatch.setattr(reference, "_download", lambda key: (None, "no reference available (offline)"))
    result = render(RenderRequest(smiles=CAMPHOR, reference=True))
    assert any("polycyclic" in note for note in result.notes)


def test_reference_off_by_default(cache):
    result = render(RenderRequest(smiles=CAMPHOR))
    assert result.depiction.startswith("automatic")


def test_reference_does_not_break_verification():
    """The reference path changes coordinates only, so the gate still applies."""
    result = render(RenderRequest(smiles=CAMPHOR, reference=True))
    assert result.report.ok
    assert "formula, charge and radical count survive preparation" in result.report.checks_run


