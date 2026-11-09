# Validation

How much of the tool's claim survives contact with real data, measured rather
than asserted.

Run: `python tools/validate.py --identity 4000 --render 700 --mutation 200`
Raw output: `data/report/validation.json`, failure lists as CSV beside it.

## The claim splits in two, and only one is testable against a dataset

1. **Identity.** The molecule parsed and drawn is the molecule that was given.
   This is testable at scale, with an external oracle.
2. **Depiction.** The picture actually shows that molecule, not merely a
   molecule derived from it. No dataset can test this. It is measured by
   corrupting drawings and counting whether the gate notices.

## Method

Data: **ChEMBL**, via the public API, 4000 small molecules. Each record carries
a canonical SMILES and an InChIKey computed by InChI, a different
implementation from the one that draws. That independence is the point.

For every molecule:

- compute the InChIKey with the reference implementation, from the dataset's
  own SMILES. If that disagrees with the dataset's stored InChIKey, the record
  is excluded, not counted as a pass or a failure. Grading against an oracle
  that disagrees with itself measures nothing.
- run our parse, then our 2D preparation, and compare InChIKeys at both stages.
- conformers are stripped before computing keys, because InChI perceives
  cis/trans from 2D coordinates. Leaving them in turns a layout choice into an
  apparent change of substance.

## Results

### Identity, n = 4000

| outcome | count |
|---|---|
| identical after parse and after drawing preparation | **3982 (99.55%, CI 99.29–99.72)** |
| excluded, dataset InChIKey disagrees with InChI on its own SMILES | 18 |
| our parse changed the molecule | 0 |
| our preparation changed the molecule | 0 |

All 18 exclusions are macrocyclic quaternary pyridinium dyes where ring
double-bond stereo is perceived differently by InChI builds. The oracle is
inconsistent there, so those records cannot judge anything.

### Gate detection, 200 valid drawings per class

| corruption | caught | rate |
|---|---|---|
| bond removed from the drawing | 200/200 | 100% |
| atom label erased | 200/200 | 100% |
| canvas too small (clipping) | 200/200 | 100% |
| charge changed after the check | 200/200 | 100% |
| wedge removed from a stereocentre | 42/42 | 100% |
| E/Z flipped | 12/12 | 100% |
| clean drawings wrongly rejected | 0/200 | 0% (CI 0–1.88) |

E/Z flips on double bonds whose reference substituent sits in a ring are not
counted. Reflecting a ring atom does not produce an alternative isomer, it
produces a broken ring, which no renderer would emit. Those cases are covered
by the identity test instead.

### Coverage, n = 700 real molecules

| outcome | count |
|---|---|
| rendered and passed verification | **700/700 (CI 99.45–100)** |
| refused | 0 |
| failed the gate | 0 |
| caveat: unspecified double bond drawn as one isomer | 23 |
| caveat: unspecified stereocentre drawn flat | 220 |

### Edge cases absent from ChEMBL, 42 molecules

Charged species, carbocations, carbanions, carbenes, radicals, isotopes,
hypervalent S and P, boron, silicon, tin, selenium, fused rings, sugars,
metal complexes, multi-fragment salts. All 42 rendered and verified.

## What validation found in the tool itself

Three real defects, all found by testing rather than review:

1. The parser used `MolFromSmiles(sanitize=False)` and never re-derived stereo,
   so `/` and `\` were dropped and the E/Z check compared zero against zero.
2. The E/Z check guessed substituents instead of using the ones the descriptor
   refers to, producing false alarms on 7% of substituted alkenes.
3. InChI reads cis/trans from coordinates, which made the identity oracle
   report changes that were layout, not chemistry.

## What is still not measured

The gate verifies the drawing structurally, by inspecting what was drawn. It
does not read the raster back and identify the molecule from pixels. A wrong
label glyph would be invisible to it, because SVG atom labels are vector
outlines, not text.

Closing that gap needs an optical structure recogniser on a sample, for example
DECIMER or MolScribe: render, read the PNG back, compare InChIKey. That is
independent of both RDKit and this gate, and it introduces its own error rate,
which would have to be reported alongside. It is the next test to run, not a
substitute for the ones above.

## Reproduce

```
pip install -e .
python tools/validate.py --identity 4000 --render 700 --mutation 200
python -m pytest tests -q
```

First run downloads and caches the ChEMBL sample in `data/chembl-sample.jsonl`.
