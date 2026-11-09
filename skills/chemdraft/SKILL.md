---
name: chemdraft
description: Draw chemistry figures as PNG or SVG. Use for molecules, Lewis structures with lone pairs, reaction schemes, mechanism arrows, comparison tables, energy diagrams, and Newman or Fischer projections. Every figure is verified against the molecule before it is written, and nothing is written when it cannot be verified.
---

# chemdraft

A command line tool. It draws, checks, and only then writes.

```
chemdraft smiles "CC(=O)O" -o acid.png
chemdraft render scheme.json -o scheme.png
```

Always quote the SMILES. PowerShell treats `[`, `]`, `\`, `$` and `+` as syntax.

## Using the output

JSON on stdout:

```json
{
  "ok": true,
  "formula": "C2H4O2",
  "mw": 60.05,
  "charge": 0,
  "stereocentres": [],
  "png": "acid.png",
  "depiction": "reference depiction (pubchem)",
  "embed": "![[acid.png]]",
  "verified": true,
  "warnings": []
}
```

Read `formula` and quote it. If it is not what you expected, you wrote the wrong
SMILES and the figure is wrong even though the tool is happy.

| exit | meaning |
|---|---|
| 0 | written and verified |
| 2 | the SMILES or spec is invalid, or a reaction does not balance |
| 3 | the drawing failed verification |
| 4 | internal error, nothing written |

On 2, 3 or 4 there is no file. Do not describe a figure that does not exist.

## Say what you cannot do

This tool draws 2D schematics and verifies them. Outside that, say so instead of
inventing something plausible.

- **No file means no figure.** If a call fails, report the error. Do not describe
  the picture you intended, and do not write an image yourself.
- **Never bend the input to pass a check.** If a `reaction` does not balance, the
  chemistry is wrong, not the tool. Do not delete a reagent or edit a SMILES
  until it goes through. Report the mismatch.
- **Do not guess a SMILES to get a result.** If you are unsure of a structure,
  say so. A wrong SMILES produces a confident, verified, wrong figure; the
  `formula` in the output is the only cross-check you have.
- **Do not draw stereochemistry the input did not specify.** If the tool reports
  unspecified centres, say the configuration is undefined rather than picking
  one.
- **Do not present one outcome as the only one.** Racemates and mixtures get
  captioned as such.
- **`caveats` are part of the answer.** If the tool says the depiction is
  automatic on a bridged skeleton, or that items overlap, pass that on instead
  of dropping it.
- **Out of scope, say so:** 3D geometry, orbitals, crystal lattices, spectra,
  and any figure whose chemistry you cannot state as SMILES. Use the primitives
  for schematic shapes if that is honest, and name the abstraction.

Being told "the tool cannot do this" is a useful answer. A figure that looks
right and is not is the one failure this tool exists to prevent.

## One molecule

```
chemdraft smiles "<smiles>" -o out.png
      [--mode skeletal|lewis] [--highlight "<SMARTS>"] [--svg] [--color]
      [--no-reference] [--bond-px N]
```

- `--mode lewis` draws every hydrogen, labels every atom, and adds lone pairs.
  Use it for bonding and electron questions. Use skeletal for anything larger
  than a few atoms.
- `--highlight` takes SMARTS and marks the match, which is how you point at a
  reactive site instead of describing it.
- `--no-reference` forces automatic layout instead of a published depiction.

## Scenes

```
chemdraft render spec.json -o out.png [--svg] [--no-reference]
```

```json
{
  "canvas": { "padding": 24, "background": "#ffffff" },
  "style":  { "font_size": 16, "bond_px": 40, "color": false, "reference": true },
  "items": [ ... ]
}
```

`items` is a list laid out top to bottom. Every item has a `type`.

### Containers

| type | fields | meaning |
|---|---|---|
| `row` | `items`, `gap`, `align` (top/center/bottom) | horizontal |
| `col` | `items`, `gap`, `align` (left/center/right) | vertical |
| `grid` | `items`, `cols`, `gap` | uniform cells |
| `panel` | `items`, `title`, `padding`, `layout` (row/col) | bordered group |
| `table` | `rows` (list of lists of items), `gap`, `border` | cells may be text or molecules |
| `reaction` | `from`, `to`, `reagents`, `above`, `below`, `length` | see balance, below |

### Leaves

| type | fields |
|---|---|
| `molecule` | `smiles`, `mode`, `caption`, `highlight`, `bond_px`, `id` |
| `text` | `content`, `size`, `align` (start/center), `italic` |
| `arrow` | `above`, `below`, `length` |
| `curved_arrow` | `from`, `to`, `curve`, `standoff`, `label` |
| `newman` | `smiles`, `bond`, `dihedral`, `bond_px`, `font_size` |
| `fischer` | `smiles`, `bond_px`, `font_size` |
| `line` | `w`, `h`, `x1`, `y1`, `x2`, `y2`, `arrow` |
| `polyline` | `w`, `h`, `points`, `smooth`, `stroke_width`, `dash` |
| `polygon` | `w`, `h`, `points`, `fill` |
| `circle` | `w`, `h`, `r`, `fill` |
| `rect` | `w`, `h`, `fill`, `radius`, `stroke`, `dash` |

Primitives draw in local coordinates from `(0,0)` at the top left of their box.

## Reaction balance

`reaction` is the correct way to draw a step, because it refuses one that does
not conserve atoms:

```json
{ "type": "reaction", "above": "HBr", "reagents": ["Br"],
  "from": ["CC=C1CCC1"], "to": ["CC1(Br)CCCC1"] }
```

`from` and `to` are lists of SMILES, drawn with `+` between them. `reagents`
count toward the left side: `HBr` is `["Br"]`, a proton is `["[H+]"]`, water is
`["O"]`. A mismatch is refused with the elements named:

```
reaction does not balance by element count: C: 6 on the left, 7 on the right
```

Use it whenever a step actually conserves atoms. Use `row` with an `arrow` for
mechanism steps, where the species are ions and a balance would be meaningless.

## Mechanism arrows

Give a molecule an `id`, then point at an atom or a bond by index. Indices are
RDKit's, in SMILES order, so `C[C@@H](N)C(=O)O` numbers the methyl 0, the
stereocentre 1, N 2, the carboxyl carbon 3.

```json
{ "id": "a", "type": "molecule", "smiles": "CC(=O)[O-]" },
{ "type": "curved_arrow",
  "from": { "item": "a", "atom": 3 },
  "to":   { "item": "a", "atom": 1 },
  "curve": 0.45 }
```

`curve` bends the arrow, negative bends the other way. `standoff` keeps the head
clear of a label. An arrow to a bond points at its midpoint: `{ "bond": 1 }`.

## Projections

```json
{ "type": "newman", "smiles": "CCCC", "bond": [1, 2], "dihedral": 60 }
{ "type": "fischer", "smiles": "C[C@H](N)C(=O)O" }
```

`newman` views along a bond; omit `bond` and it takes the most substituted C–C
bond. `dihedral` is 180 anti, 60 gauche, 0 eclipsed.

`fischer` needs an open chain of at least three carbons and puts the most
oxidised end at the top. Choose the enantiomer through the SMILES: `C[C@H](N)C(=O)O`
is L-alanine, `C[C@@H](N)C(=O)O` is D.

## When a figure needs chemistry the tool cannot know

The gate checks the drawing against the SMILES. It cannot know what the flask
contains. State that yourself:

- **Not enantioselective** (hydroboration, most additions): caption the product
  `"(1R,2R)-2-methylcyclohexanol (± racemate)"`.
- **One outcome of several**: write "major" or name the alternative.
- **Unspecified stereocentres**: the tool reports them as drawn flat; say so if
  it matters.
- **A product that disagrees with its starting material by an atom**: state both
  as SMILES in a `reaction` item and let the tool check, or the scheme will look
  convincing and be impossible.

## Recipes

**Acid with lone pairs**
```json
{ "type": "molecule", "smiles": "CC(=O)O", "mode": "lewis", "caption": "acetic acid" }
```

**Two-step scheme**
```json
{ "type": "col", "gap": 24, "items": [
  { "type": "reaction", "above": "H\u2082SO\u2084, \u0394",
    "from": ["CC(=O)O", "CCO"], "to": ["CC(=O)OCC", "O"] },
  { "type": "reaction", "above": "NaOH", "below": "hydrolysis",
    "from": ["CC(=O)OCC", "O"], "to": ["CC(=O)O", "CCO"] }
]}
```

Organometallic and other steps whose byproducts are usually ignored will not
balance on their own. Either include the byproducts (`Mg(OH)Br` is `"O[Mg]Br"`)
or draw the step as a `row` with an `arrow`, where no balance is claimed.

**Comparison table** (name / formula / Lewis / line)
```json
{ "type": "table", "border": true, "gap": 20, "rows": [
  [{ "type": "text", "content": "NAME" }, { "type": "text", "content": "LINE" }],
  [{ "type": "text", "content": "ethanol" }, { "type": "molecule", "smiles": "CCO" }]
]}
```

**Energy diagram**
```json
{ "type": "col", "gap": 6, "items": [
  { "type": "polyline", "w": 420, "h": 150, "smooth": true,
    "points": [[10,130],[70,40],[130,95],[280,30],[410,120]] },
  { "type": "line", "w": 420, "h": 20, "x1": 10, "y1": 10, "x2": 410, "y2": 10,
    "arrow": true },
  { "type": "text", "content": "reaction coordinate", "size": 14 }
]}
```

## What is verified

- the SMILES parses and is chemically possible
- formula, charge and radical count are unchanged from input to drawing
- every bond appears; every atom that needs a symbol has one
- assigned stereocentres carry a wedge that **encodes** the stated configuration
- double bonds are drawn cis or trans as declared
- no atom is drawn outside the canvas or on top of another
- in a scene, every molecule passes the same checks, and `reaction` items balance

Overlapping items are reported in `caveats`, not as failures.

## Limits

- Faithful to the SMILES, not to intent. A wrong SMILES yields a confident wrong
  figure; the formula in the output is the cross-check.
- Automatic layout is chosen by geometric score. For a bridged skeleton without
  a published depiction it can be unfamiliar; `depiction` says when that is the
  case.
- Fischer takes the longest carbon chain and refuses rings. Newman needs a bond
  with substituents at both ends.
- No 3D geometry, orbitals or crystal structures. Use primitives for schematic
  diagrams of those.

## Putting figures in a course vault

```
chemdraft smiles "CCO" -o "<vault>/Courses/Chemistry/visuals/ethanol.png"
```

Then embed the `embed` value from the output in the topic note. Keep figures as
files under `visuals/`, not pasted into chat.
