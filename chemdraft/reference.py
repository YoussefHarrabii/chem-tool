"""Reference depictions.

Automatic layout optimises geometry. It has no idea how a chemist draws a
given skeleton, and on bridged systems it produces valid but unfamiliar
projections, which is a problem when the figure is meant to teach someone who
will meet the conventional drawing in an exam.

So: ask PubChem for a published 2D depiction of this exact molecule, by
InChIKey. Verify that what came back really is the same molecule, then use its
coordinates. Offline, unknown, or mismatched: fall back to automatic layout and
say so, instead of pretending the result is canonical.

A record is about 4 KB. Cached on disk, so a molecule is fetched once.
"""

from __future__ import annotations

import os
import urllib.error
import urllib.request
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import inchi

CACHE = Path(os.environ.get("LOCALAPPDATA") or Path.home() / ".cache") / "chemdraft" / "refs"
TIMEOUT = 8
API = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{key}"
       "/record/SDF?record_type=2d")


def cache_path(key: str) -> Path:
    return CACHE / f"{key}.sdf"


def _download(key: str) -> tuple[str | None, str]:
    url = API.format(key=key)
    try:
        with urllib.request.urlopen(url, timeout=TIMEOUT) as response:
            return response.read().decode("utf-8", "replace"), "pubchem"
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            return None, "not in PubChem"
        return None, f"PubChem returned HTTP {exc.code}"
    except Exception as exc:  # noqa: BLE001 - network trouble must not break drawing
        return None, f"no reference available ({type(exc).__name__})"


def record_for(mol: Chem.Mol) -> tuple[Chem.Mol | None, str]:
    """A verified reference molecule with 2D coordinates, or (None, reason)."""
    try:
        key = inchi.MolToInchiKey(mol)
    except Exception:  # noqa: BLE001
        return None, "no InChIKey"
    if not key:
        return None, "no InChIKey"

    path = cache_path(key)
    if path.exists():
        block, source = path.read_text(encoding="utf-8", errors="replace"), "cache"
    else:
        block, source = _download(key)
        if block is None:
            return None, source
        try:
            CACHE.mkdir(parents=True, exist_ok=True)
            path.write_text(block, encoding="utf-8")
        except OSError:
            pass  # a read-only cache is not a reason to fail

    ref = Chem.MolFromMolBlock(block, removeHs=True, sanitize=True)
    if ref is None:
        return None, "reference record could not be read"

    # The record must be the same molecule. PubChem is authoritative for what it
    # holds, but a wrong or stale record must never be drawn.
    try:
        if inchi.MolToInchiKey(ref) != key:
            return None, "reference record does not match this molecule"
    except Exception:  # noqa: BLE001
        return None, "reference record has no usable InChIKey"
    return ref, source


def apply_coordinates(mol: Chem.Mol, ref: Chem.Mol) -> bool:
    """Copy the reference coordinates onto mol. False when the atom sets do not
    correspond, in which case the caller keeps its own layout."""
    if ref.GetNumAtoms() != mol.GetNumAtoms():
        return False
    if ref.GetNumConformers() == 0:
        return False
    match = mol.GetSubstructMatch(ref)
    if len(match) != ref.GetNumAtoms():
        return False

    ref_conf = ref.GetConformer()
    conf = Chem.Conformer(mol.GetNumAtoms())
    for ref_idx, mol_idx in enumerate(match):
        p = ref_conf.GetAtomPosition(ref_idx)
        conf.SetAtomPosition(mol_idx, p)
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    return True


def polycyclic(mol: Chem.Mol) -> bool:
    """True when an atom belongs to more than one ring.

    Polycyclic skeletons are where automatic layout most often picks a
    projection a chemist would not, because the flat drawing is a choice among
    several and geometry alone does not pick the familiar one.
    """
    info = mol.GetRingInfo()
    return any(info.NumAtomRings(i) > 1 for i in range(mol.GetNumAtoms()))
