"""Parse and describe molecules. Never repairs, never guesses.

The contract: if this module returns a molecule, that molecule is exactly what
the SMILES said. If the SMILES is invalid or chemically impossible, this module
raises and nothing downstream runs.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors

# RDKit prints parse complaints to stderr. We surface them as exceptions
# instead, so a caller can never mistake a partial parse for a success.
RDLogger.DisableLog("rdApp.*")


class ChemError(Exception):
    """The input is not a molecule we can draw faithfully."""


@dataclass
class Description:
    smiles: str
    canonical: str
    formula: str
    mw: float
    charge: int
    atoms: int
    bonds: int
    stereocentres: list[str] = field(default_factory=list)
    unassigned_stereocentres: int = 0
    double_bond_stereo: int = 0
    radicals: int = 0


def parse(smiles: str) -> Chem.Mol:
    """SMILES to Mol. Raises ChemError with a readable reason, never returns
    None and never falls back to a partially sanitised structure."""
    if not smiles or not smiles.strip():
        raise ChemError("empty SMILES")

    raw = Chem.MolFromSmiles(smiles, sanitize=False)
    if raw is None:
        raise ChemError(
            f"not valid SMILES: {smiles!r}. Check ring closures, bracket "
            f"balance and atom valences."
        )

    # DetectChemistryProblems reports valence and aromaticity faults with
    # atom indices before sanitisation sanitises them away.
    problems = Chem.DetectChemistryProblems(raw)
    if problems:
        msgs = []
        for p in problems:
            idx = p.GetAtomIdx() if hasattr(p, "GetAtomIdx") else None
            where = f"atom {idx}" if idx is not None else "molecule"
            msgs.append(f"{where}: {p.Message()}")
        raise ChemError("chemically impossible structure, " + "; ".join(msgs))

    try:
        Chem.SanitizeMol(raw)
    except Exception as exc:  # noqa: BLE001 - we re-raise with context
        raise ChemError(f"structure failed sanitisation: {exc}") from exc

    # Sanitisation alone does not derive stereochemistry. Without these two
    # calls, stereo written as @, @@, / or \ is silently dropped and the figure
    # would claim a stereochemistry the input never had.
    Chem.AssignStereochemistry(raw, cleanIt=True, force=True)
    Chem.SetBondStereoFromDirections(raw)
    Chem.AssignStereochemistry(raw, cleanIt=True, force=True)

    # Round-trip invariant: canonicalise, reparse, canonicalise again. If the
    # two disagree, something in the parse lost information and we refuse.
    canonical = Chem.MolToSmiles(raw)
    again = Chem.MolFromSmiles(canonical)
    if again is None or Chem.MolToSmiles(again) != canonical:
        raise ChemError(
            "this structure does not survive a canonical round trip, so it "
            "cannot be drawn faithfully"
        )

    return raw


def describe(mol: Chem.Mol, smiles: str = "") -> Description:
    centres = Chem.FindMolChiralCenters(
        mol, includeUnassigned=True, useLegacyImplementation=False
    )
    assigned = [f"{idx}:{cfg}" for idx, cfg in centres if cfg != "?"]
    unassigned = sum(1 for _, cfg in centres if cfg == "?")
    db_stereo = sum(
        1 for b in mol.GetBonds() if b.GetStereo() != Chem.BondStereo.STEREONONE
    )
    radicals = sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
    return Description(
        smiles=smiles,
        canonical=Chem.MolToSmiles(mol),
        formula=rdMolDescriptors.CalcMolFormula(mol),
        mw=round(Descriptors.MolWt(mol), 2),
        charge=Chem.GetFormalCharge(mol),
        atoms=mol.GetNumAtoms(),
        bonds=mol.GetNumBonds(),
        stereocentres=assigned,
        unassigned_stereocentres=unassigned,
        double_bond_stereo=db_stereo,
        radicals=radicals,
    )


_PT = Chem.GetPeriodicTable()


def bond_order_sum(atom: Chem.Atom) -> float:
    total = 0.0
    for bond in atom.GetBonds():
        kind = bond.GetBondTypeAsDouble()
        if kind == 1.5:  # aromatic, kept as 1.5 (a kekulised copy has integers)
            total += 1.5
        else:
            total += kind
    return total


def needs_label(atom: Chem.Atom) -> bool:
    """True when the drawing must show a symbol, or it would describe a
    different species.

    A bare vertex is honest only when the viewer infers exactly this atom from
    it: carbon, uncharged, with the hydrogen count that filling its valence
    would produce. Everything else must be written out.
    """
    if atom.GetSymbol() != "C":
        return True
    if atom.GetFormalCharge() != 0:
        return True
    if atom.GetIsotope() != 0:
        return True
    if atom.GetNumRadicalElectrons() != 0:
        return True
    valence = _PT.GetDefaultValence(atom.GetAtomicNum())
    implied_h = max(0.0, valence - bond_order_sum(atom))
    if abs(atom.GetTotalNumHs() - implied_h) > 1e-6:
        return True
    return False


def kekulized_copy(mol: Chem.Mol) -> Chem.Mol:
    """Aromatic bonds of order 1.5 break valence arithmetic, so count electrons
    on an integer-bond copy."""
    copy = Chem.Mol(mol)
    try:
        Chem.Kekulize(copy, clearAromaticFlags=True)
    except Exception:  # noqa: BLE001 - keep aromatic if kekulisation fails
        return Chem.Mol(mol)
    return copy
