# chemdraft

Draw chemistry figures from the command line, checked before they are written.

```
chemdraft smiles "CC(=O)O" -o acid.png
```

Molecules, Lewis structures with lone pairs, reaction schemes, mechanism arrows,
tables, energy diagrams, Newman and Fischer projections. Input is SMILES or a
small JSON scene; output is PNG and SVG.

Built for a coding agent working with a student: the figure is verified against
the molecule before it reaches disk, and when it cannot be verified no file is
written at all.

## Install

```
pip install -e .
```

Python 3.10+, RDKit, PyMuPDF, Pillow. No TeX, no headless browser, no server.

## Two entry points

**One molecule.** `chemdraft smiles "<smiles>" -o out.png` with optional
`--mode lewis`, `--highlight "<SMARTS>"`, `--svg`, `--color`, `--no-reference`.

**A scene.** `chemdraft render spec.json -o out.png` for reactions, schemes,
tables, projections. The format is documented in `skills/chemdraft/SKILL.md`,
and `examples/` has a working spec for each feature.

Every call returns JSON on stdout with the formula, molecular weight, charge,
stereocentres, the depiction used, warnings, and the paste-ready
`![[name.png]]` embed string.

## What is checked before a file exists

The molecule is parsed without repair. Validity, valence, formula, charge and
stereochemistry are compared between the input and the molecule that is drawn.
Every bond must appear, every atom that needs a symbol must have one, every
assigned stereocentre must carry a wedge that encodes its configuration, and
double bonds must be drawn cis or trans as declared. In a scene, an arrow
between `from` and `to` must conserve elements.

On failure: exit code 2 for bad input, 3 for a failed drawing, 4 for an internal
error. No file is written in any of those cases.

## Reference depictions

Automatic layout optimises geometry but has no idea how a chemist draws a given
skeleton, and on bridged systems it produces valid but unfamiliar projections.
So the tool asks PubChem for a published 2D depiction by InChIKey, verifies the
record is the same molecule, and uses its coordinates. A record is about 4 KB
and is cached under `%LOCALAPPDATA%\chemdraft\refs`. Offline or unknown, it
falls back to automatic layout and says so.

## Layout

```
chemdraft/
  chemistry.py   parse, describe, sanitise, lone pair arithmetic
  render.py      single molecule: canvas, draw, verify, rasterise
  verify.py      the gate
  layout.py      choose the best of several 2D layouts by geometric score
  annotate.py    lone pair dots
  reference.py   PubChem depictions, cached and verified
  projection.py  Newman and Fischer, verified by reading the convention back
  scene.py       spec, layout, composition, atom balance, overlap lint
  cli.py         smiles / render
tests/           70 tests, including deliberately broken drawings
tools/           validate.py, verify_sheet.py, quiz_sheet.py
examples/        one working spec per feature
```

## Design rules

- **No special cases.** No branch in the production code keys on a particular
  molecule, figure, or element's display. Where something must be known it comes
  from chemistry data or from the spec.
- **Fail closed.** A missing figure is safe; a wrong one teaches the wrong thing.
- **Verify, do not assume.** Stereochemistry is read back from the drawing, the
  Fischer arrangement is re-derived from its own geometry, and Newman angles are
  asserted numerically rather than eyeballed.

## Validation

`VALIDATION.md` records the measurements: 4,000 ChEMBL molecules for identity
against a published InChIKey (99.55%, zero failures of our own), 1,200
deliberately corrupted drawings for gate detection, and the known gaps.

## Known limits

- Faithful to the SMILES, not to your intent. A wrong SMILES produces a
  confident wrong figure; the formula in the output is the cross-check.
- Reaction layout is not automatic. A long scheme is placed by its spec.
- Fischer takes the longest carbon chain and refuses rings. Newman needs a bond
  with substituents on both ends.
- No 3D geometry, orbitals or crystal structures.
