"""Measure the tool instead of trusting it.

Three separate questions, because they fail for different reasons:

1. Identity  - does our pipeline keep the molecule the dataset says it is?
               Oracle: ChEMBL's InChIKey, computed by InChI, not by us.
2. Drawing   - does the gate catch a drawing that no longer shows the molecule?
               Measured by corrupting good drawings on purpose and counting
               detections. A gate nobody tries to break is worth nothing.
3. Coverage  - what fraction of real molecules render and pass? Failures are
               written out so they can be classified rather than averaged away.

    python tools/validate.py --identity 10000 --render 1200

Results go to data/report/.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import urllib.request
from collections import Counter
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import inchi as rdinchi

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from chemdraft.chemistry import ChemError, kekulized_copy, parse  # noqa: E402
from chemdraft.render import RenderRequest, render  # noqa: E402
from chemdraft.verify import verify  # noqa: E402

DATA = ROOT / "data"
REPORT = DATA / "report"
CACHE = DATA / "chembl-sample.jsonl"
API = ("https://www.ebi.ac.uk/chembl/api/data/molecule.json"
       "?limit=1000&offset={offset}")

# Molecules ChEMBL does not contain: charged intermediates, radicals, isotopes,
# hypervalent centres, main-group elements. The distribution medicinal
# chemistry never touches but a chemistry lesson does.
EDGE_CASES = [
    "[CH3+]", "[CH3-]", "[CH2]", "C[CH+]C", "C[CH-]C",
    "[C-]#[O+]", "[O-][N+](=O)C", "C[N+](C)(C)C", "[NH4+]",
    "[13CH4]", "[2H]O[2H]", "C[13C](=O)O",
    "O=S(=O)(O)O", "CS(=O)C", "CP(=O)(O)O", "O=[P-]([O-])[O-]",
    "OB(O)c1ccccc1", "C[Si](C)(C)C", "C[Sn](C)(C)C", "[SeH2]",
    "c1cc2ccc3cccc4ccc(c1)c2c34", "C1CC2CCC1CC2",
    "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
    "C[C@@H](O)[C@H](N)C(=O)O", "F/C=C/F", "F/C=C\\F",
    "C/C=C/C=C/C", "N#[N+][O-]", "c1ccc([nH]1)C", "c1ccoc1", "c1ccsc1",
    "[Na+].[Cl-]", "CC(=O)[O-].[Na+]",
    "C[C@H]1CC[C@@H](C)CC1", "O=[Mn](=O)(=O)[O-]", "[Fe+2]",
    "C(#N)C#N", "N=C=O", "O=C=O", "ClC(Cl)(Cl)Cl", "FC(F)(F)F",
    "C12C3C4C1C5C2C3C45",
]


# ---------------------------------------------------------------- fetching ---
def load_sample(n: int) -> list[dict]:
    DATA.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    if CACHE.exists():
        with CACHE.open(encoding="utf-8") as fh:
            for line in fh:
                rows.append(json.loads(line))
        if len(rows) >= n:
            return rows[:n]

    offset = 0
    t0 = time.time()
    while len(rows) < n:
        url = API.format(offset=offset)
        with urllib.request.urlopen(url, timeout=120) as resp:
            payload = json.loads(resp.read().decode("utf-8"))
        batch = payload.get("molecules") or []
        if not batch:
            break
        for m in batch:
            struct = m.get("molecule_structures") or {}
            smiles = struct.get("canonical_smiles")
            key = struct.get("standard_inchi_key")
            if not smiles or not key:
                continue
            props = m.get("molecule_properties") or {}
            rows.append({
                "id": m.get("molecule_chembl_id"),
                "smiles": smiles,
                "ref_key": key,
                "ref_formula": props.get("full_molformula"),
                "ref_mwt": props.get("full_mwt"),
            })
        offset += 1000
        print(f"  fetched {len(rows)}  ({time.time() - t0:.0f}s)", flush=True)

    with CACHE.open("w", encoding="utf-8") as fh:
        for r in rows[:n]:
            fh.write(json.dumps(r) + "\n")
    return rows[:n]


# ------------------------------------------------------------- primitives ---
def key_of(mol: Chem.Mol | None) -> str | None:
    """InChIKey with conformers stripped.

    InChI perceives cis/trans from 2D coordinates, so leaving the drawn
    coordinates in the molecule turns a layout decision into an apparent change
    of substance. Strip them, and the key measures connectivity and declared
    stereochemistry only, which is what we are testing.
    """
    if mol is None:
        return None
    try:
        flat = Chem.Mol(mol)
        flat.RemoveAllConformers()
        k = rdinchi.MolToInchiKey(flat)
        return k or None
    except Exception:  # noqa: BLE001
        return None


def prepared_of(mol: Chem.Mol, lewis: bool = False) -> Chem.Mol:
    work = Chem.Mol(mol)
    if lewis:
        work = Chem.AddHs(work)
    rdDepictor.Compute2DCoords(work)
    return rdMolDraw2D.PrepareMolForDrawing(work)


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, centre - half), min(1.0, centre + half))


def rate(name: str, k: int, n: int) -> str:
    lo, hi = wilson(k, n)
    return f"{name:38s} {k:6d}/{n:<6d} {100 * k / max(n, 1):6.2f}%  CI {100 * lo:.2f}-{100 * hi:.2f}%"


# ------------------------------------------------------------- test 1 ------
def test_identity(rows: list[dict]) -> dict:
    print("\n[1] identity, oracle = ChEMBL standard_inchi_key")
    stats = Counter()
    failures = []
    for r in rows:
        stats["total"] += 1
        # Baseline: the reference implementation, so a dataset artefact is not
        # blamed on our code.
        base = Chem.MolFromSmiles(r["smiles"])
        base_key = key_of(base)
        if base_key is None:
            stats["oracle_cannot_parse"] += 1
            continue
        if base_key != r["ref_key"]:
            stats["dataset_key_disagrees"] += 1
            continue
        try:
            ours = parse(r["smiles"])
        except ChemError as exc:
            stats["we_refused"] += 1
            failures.append({**r, "stage": "parse", "detail": str(exc)})
            continue
        ours_key = key_of(ours)
        if ours_key != base_key:
            stats["parse_changed_molecule"] += 1
            failures.append({**r, "stage": "parse_key", "detail": f"{ours_key} != {base_key}"})
            continue
        try:
            prep = prepared_of(ours)
        except Exception as exc:  # noqa: BLE001
            stats["prepare_crashed"] += 1
            failures.append({**r, "stage": "prepare", "detail": repr(exc)})
            continue
        prep_key = key_of(prep)
        if prep_key != base_key:
            stats["prepare_changed_molecule"] += 1
            failures.append({**r, "stage": "prepare_key", "detail": f"{prep_key} != {base_key}"})
            continue
        stats["exact"] += 1

    n = stats["total"]
    print(rate("identity preserved end to end", stats["exact"], n))
    for k, v in sorted(stats.items()):
        if k not in ("total", "exact") and v:
            print(f"    {k:36s} {v}")
    with (REPORT / "identity-failures.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["id", "smiles", "stage", "detail", "ref_key"], extrasaction="ignore")
        w.writeheader()
        w.writerows(failures[:500])
    return {"stats": dict(stats), "failures": len(failures)}


# ------------------------------------------------------------- test 2 ------
def mutate_cases(result) -> tuple[list[tuple[str, object]], list[str]]:
    """Corrupt a good drawing in ways that change what it says, and list the
    lies the gate is expected to catch."""
    out: list[tuple[str, object]] = []
    skipped: list[str] = []
    prepared = result.prepared

    def recheck(**kw):
        args = dict(svg=result.svg, width=result.width, height=result.height,
                    label_all=False)
        args.update(kw)
        return verify(parsed=prepared, prepared=prepared, drawer=result.drawer,
                      expected=result.description, **args)

    svg = result.svg
    dropped = svg
    for b in range(prepared.GetNumBonds()):
        dropped = dropped.replace(f"class='bond-{b}", f"class='gone-{b}")
    out.append(("bond removed from drawing", lambda: recheck(svg=dropped)))

    stripped = svg
    for a in prepared.GetAtoms():
        stripped = stripped.replace(f"class='atom-{a.GetIdx()}'", f"class='gone-{a.GetIdx()}'")
    out.append(("label erased from drawing", lambda: recheck(svg=stripped)))

    def clip():
        return recheck(width=12, height=12)
    out.append(("canvas too small", clip))

    if any(b.GetStereo() != Chem.BondStereo.STEREONONE for b in prepared.GetBonds()):
        stereo_bond = next(b for b in prepared.GetBonds()
                           if b.GetStereo() != Chem.BondStereo.STEREONONE)
        refs = list(stereo_bond.GetStereoAtoms())
        rings = prepared.GetRingInfo()
        if refs and any(rings.NumAtomRings(a) > 0 for a in refs):
            # Reflecting a ring atom does not produce an alternative drawing,
            # it produces a broken ring. No renderer would do that, so this
            # corruption is not a fair test and is skipped rather than scored.
            skipped.append("E/Z flipped (reference atom is in a ring)")
        else:
            def flip():
                conf = prepared.GetConformer()
                bond = [b for b in prepared.GetBonds()
                        if b.GetStereo() != Chem.BondStereo.STEREONONE][0]
                i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                nbr = [n.GetIdx() for n in prepared.GetAtomWithIdx(j).GetNeighbors()
                       if n.GetIdx() != i]
                if not nbr:
                    return recheck()
                k = nbr[0]
                pi, pj, p = (np.array([conf.GetAtomPosition(x).x,
                                       conf.GetAtomPosition(x).y])
                             for x in (i, j, k))
                u = (pj - pi) / np.linalg.norm(pj - pi)
                v = p - pi
                perp = v - np.dot(v, u) * u
                if np.linalg.norm(perp) < 1e-3:
                    return recheck()
                f = pi + v - 2 * perp
                conf.SetAtomPosition(
                    k, Chem.rdGeometry.Point3D(float(f[0]), float(f[1]), 0.0)
                )
                return recheck()
            out.append(("E/Z flipped", flip))

    if any(cfg != "?" for _, cfg in Chem.FindMolChiralCenters(
            prepared, includeUnassigned=True, useLegacyImplementation=False)):
        def unwind():
            for bond in prepared.GetBonds():
                bond.SetBondDir(Chem.BondDir.NONE)
            return recheck()
        out.append(("wedge removed from stereocentre", unwind))

    out.append(("charge silently changed", lambda: _charge_mutated(prepared, recheck)))
    return out, skipped


def _charge_mutated(prepared, recheck):
    atom = next(a for a in prepared.GetAtoms() if a.GetSymbol() == "C")
    atom.SetFormalCharge(1)
    try:
        return recheck()
    finally:
        atom.SetFormalCharge(0)


def test_gate(rows: list[dict], sample: int) -> dict:
    print("\n[2] gate detection, by corrupting good drawings on purpose")
    detected = Counter()
    attempted = Counter()
    misses = []
    false_alarms = 0
    used = 0
    for r in rows:
        if used >= sample:
            break
        try:
            result = render(RenderRequest(smiles=r["smiles"], rasterize=False))
        except ChemError:
            continue
        if not result.report.ok:
            continue
        used += 1
        if not result.report.ok:
            false_alarms += 1
        for name, apply in mutate_cases(result)[0]:
            attempted[name] += 1
            try:
                report = apply()
            except Exception:  # noqa: BLE001
                detected[name] += 1  # a crash still blocks the file
                continue
            if not report.ok:
                detected[name] += 1
            else:
                misses.append({"smiles": r["smiles"], "mutation": name})
        # restore the drawing for the next mutation class
        result = None

    print(rate("clean drawings the gate wrongly rejects", false_alarms, used))
    for name in attempted:
        print(rate(f"  caught: {name}", detected[name], attempted[name]))
    if misses:
        with (REPORT / "gate-misses.csv").open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=["smiles", "mutation"])
            w.writeheader()
            w.writerows(misses)
        print(f"    {len(misses)} misses written to gate-misses.csv")
    return {"used": used, "detected": dict(detected), "attempted": dict(attempted),
            "false_alarms": false_alarms, "misses": misses}


# ------------------------------------------------------------- test 3 ------
def test_coverage(rows: list[dict], sample: int) -> dict:
    print("\n[3] coverage, how much real data renders and passes")
    stats = Counter()
    reasons = Counter()
    failures = []
    for r in rows[:sample]:
        stats["total"] += 1
        try:
            result = render(RenderRequest(smiles=r["smiles"], rasterize=False))
        except ChemError as exc:
            stats["refused"] += 1
            reasons[str(exc).split(":")[0][:60]] += 1
            failures.append({**r, "stage": "refused", "detail": str(exc)})
            continue
        except Exception as exc:  # noqa: BLE001
            stats["crashed"] += 1
            failures.append({**r, "stage": "crash", "detail": repr(exc)})
            continue
        if not result.report.ok:
            stats["failed_gate"] += 1
            first = result.report.problems[0][:80]
            reasons[first] += 1
            failures.append({**r, "stage": "gate", "detail": "; ".join(result.report.problems)})
            continue
        stats["rendered"] += 1
        if any("no specified geometry" in n for n in result.report.notes):
            stats["caveat: unspecified double bond"] += 1
        if any("drawn flat" in n for n in result.report.notes):
            stats["caveat: unspecified stereocentre"] += 1
        if result.description.formula != r["ref_formula"]:
            stats["formula_differs"] += 1
            failures.append({**r, "stage": "formula",
                             "detail": f"{result.description.formula} vs {r['ref_formula']}"})

    n = stats["total"]
    print(rate("rendered and verified", stats["rendered"], n))
    for k, v in sorted(stats.items()):
        if k not in ("total", "rendered") and v:
            print(f"    {k:36s} {v}")
    for reason, count in reasons.most_common(8):
        print(f"      {count:5d}  {reason}")
    with (REPORT / "coverage-failures.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["id", "smiles", "stage", "detail", "ref_formula"], extrasaction="ignore")
        w.writeheader()
        w.writerows(failures[:800])
    return {"stats": dict(stats), "reasons": dict(reasons)}


# ------------------------------------------------------------- test 4 ------
def test_edge_cases() -> dict:
    print("\n[4] edge cases the dataset does not contain")
    stats = Counter()
    for smiles in EDGE_CASES:
        stats["total"] += 1
        try:
            result = render(RenderRequest(smiles=smiles, rasterize=False))
        except ChemError as exc:
            stats["refused"] += 1
            print(f"    refused {smiles:30s} {str(exc)[:70]}")
            continue
        except Exception as exc:  # noqa: BLE001
            stats["crashed"] += 1
            print(f"    crash   {smiles:30s} {exc!r}"[:120])
            continue
        if not result.report.ok:
            stats["failed_gate"] += 1
            print(f"    gate    {smiles:30s} {result.report.problems[0][:70]}")
            continue
        stats["ok"] += 1
    print(rate("edge cases rendered and verified", stats["ok"], stats["total"]))
    return {"stats": dict(stats)}


# ---------------------------------------------------------------- main -----
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--identity", type=int, default=10000)
    ap.add_argument("--render", type=int, default=1200)
    ap.add_argument("--mutation", type=int, default=250)
    args = ap.parse_args()

    REPORT.mkdir(parents=True, exist_ok=True)
    print(f"loading ChEMBL sample (target {args.identity})")
    rows = load_sample(args.identity)
    print(f"  {len(rows)} molecules with SMILES and InChIKey")

    report = {
        "generated": time.strftime("%Y-%m-%d %H:%M"),
        "sample_size": len(rows),
        "identity": test_identity(rows),
        "gate": test_gate(rows, args.mutation),
        "coverage": test_coverage(rows, args.render),
        "edge": test_edge_cases(),
    }
    (REPORT / "validation.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"\nwritten to {REPORT / 'validation.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

