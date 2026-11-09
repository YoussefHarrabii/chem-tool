"""Scene composition.

A scheme must not be able to contain a structure that a single figure would
have refused, so these tests check that failure inside a scene fails the whole
scene rather than drawing a wrong panel.
"""

from __future__ import annotations

import pytest

from chemdraft.scene import SpecError, render_scene

OFFLINE = {"reference": False}


def spec(items, **style):
    return {"style": {**OFFLINE, **style}, "items": items}


def test_two_molecules_with_an_arrow():
    scene = render_scene(spec([
        {"type": "row", "items": [
            {"type": "molecule", "smiles": "CCO", "caption": "1"},
            {"type": "arrow", "above": "H2SO4", "below": "heat"},
            {"type": "molecule", "smiles": "CC=O", "caption": "2"},
        ]}
    ]))
    assert scene.molecules == 2
    assert scene.png[:8] == b"\x89PNG\r\n\x1a\n"
    assert scene.width > 0 and scene.height > 0


def test_grid_and_nesting():
    scene = render_scene(spec([
        {"type": "grid", "cols": 2, "items": [
            {"type": "molecule", "smiles": "C"},
            {"type": "molecule", "smiles": "CC"},
            {"type": "molecule", "smiles": "CCC"},
            {"type": "molecule", "smiles": "CCCC"},
        ]}
    ]))
    assert scene.molecules == 4


def test_bad_molecule_fails_the_whole_scene():
    with pytest.raises(SpecError) as exc:
        render_scene(spec([
            {"type": "row", "items": [
                {"type": "molecule", "smiles": "CCO"},
                {"type": "molecule", "smiles": "C(C)(C)(C)(C)C"},
            ]}
        ]))
    assert "failed verification" in str(exc.value) or "valence" in str(exc.value).lower()


def test_unknown_item_type_is_reported():
    with pytest.raises(SpecError) as exc:
        render_scene(spec([{"type": "blob", "smiles": "C"}]))
    assert "unknown item type" in str(exc.value)


def test_molecule_without_smiles_is_reported():
    with pytest.raises(SpecError) as exc:
        render_scene(spec([{"type": "molecule", "caption": "oops"}]))
    assert "smiles" in str(exc.value)


def test_text_item():
    scene = render_scene(spec([
        {"type": "row", "items": [
            {"type": "text", "content": "C2H6O = ethanol"},
            {"type": "molecule", "smiles": "CCO"},
        ]}
    ]))
    assert scene.molecules == 1


def test_table_with_cells():
    scene = render_scene(spec([
        {"type": "table", "border": True, "rows": [
            [{"type": "text", "content": "NAME"}, {"type": "text", "content": "LINE"}],
            [{"type": "text", "content": "ethanol"}, {"type": "molecule", "smiles": "CCO"}],
        ]}
    ]))
    assert scene.molecules == 1


def test_table_handles_ragged_rows():
    scene = render_scene(spec([
        {"type": "table", "rows": [
            [{"type": "text", "content": "a"}],
            [{"type": "text", "content": "b"}, {"type": "text", "content": "c"}],
        ]}
    ]))
    assert scene.png[:4] == b"\x89PNG"


def test_panel_with_title():
    scene = render_scene(spec([
        {"type": "panel", "title": "comparison", "items": [
            {"type": "molecule", "smiles": "CCO"},
        ]}
    ]))
    assert scene.molecules == 1


def test_primitives():
    scene = render_scene(spec([
        {"type": "col", "items": [
            {"type": "polyline", "w": 200, "h": 100, "smooth": True,
             "points": [[0, 90], [60, 20], [120, 70], [200, 30]]},
            {"type": "line", "w": 200, "h": 20, "x1": 0, "y1": 10,
             "x2": 200, "y2": 10, "arrow": True},
            {"type": "circle", "w": 60, "h": 60, "r": 25},
            {"type": "rect", "w": 80, "h": 40, "fill": "#eeeeee"},
            {"type": "polygon", "w": 60, "h": 60, "points": [[0, 60], [30, 0], [60, 60]]},
        ]}
    ]))
    assert scene.png[:4] == b"\x89PNG"


def test_primitives_need_a_size():
    with pytest.raises(SpecError) as exc:
        render_scene(spec([{"type": "line", "x1": 0, "y1": 0}]))
    assert "w" in str(exc.value)


def test_curved_arrow_between_atoms():
    scene = render_scene(spec([
        {"type": "row", "items": [
            {"id": "a", "type": "molecule", "smiles": "CC(=O)[O-]"},
            {"type": "curved_arrow",
             "from": {"item": "a", "atom": 3},
             "to": {"item": "a", "atom": 1}},
        ]}
    ]))
    assert scene.png[:4] == b"\x89PNG"


def test_curved_arrow_to_bond():
    scene = render_scene(spec([
        {"type": "row", "items": [
            {"id": "a", "type": "molecule", "smiles": "CC=O"},
            {"type": "curved_arrow",
             "from": {"item": "a", "bond": 1},
             "to": {"item": "a", "atom": 1}},
        ]}
    ]))
    assert scene.png[:4] == b"\x89PNG"


def test_curved_arrow_needs_existing_item():
    with pytest.raises(SpecError) as exc:
        render_scene(spec([
            {"type": "row", "items": [
                {"type": "molecule", "smiles": "CCO"},
                {"type": "curved_arrow",
                 "from": {"item": "missing", "atom": 0},
                 "to": {"item": "missing", "atom": 1}},
            ]}
        ]))
    assert "must name molecule items" in str(exc.value)


def test_curved_arrow_rejects_missing_atom():
    with pytest.raises(SpecError) as exc:
        render_scene(spec([
            {"type": "row", "items": [
                {"id": "a", "type": "molecule", "smiles": "CCO"},
                {"type": "curved_arrow",
                 "from": {"item": "a", "atom": 99},
                 "to": {"item": "a", "atom": 0}},
            ]}
        ]))
    assert "does not exist" in str(exc.value)


def test_duplicate_ids_are_rejected():
    with pytest.raises(SpecError) as exc:
        render_scene(spec([
            {"type": "row", "items": [
                {"id": "x", "type": "molecule", "smiles": "CCO"},
                {"id": "x", "type": "molecule", "smiles": "CCN"},
            ]}
        ]))
    assert "duplicate item id" in str(exc.value)


# ------------------------------------------------------------- reactions ---
def test_balanced_reaction_renders():
    scene = render_scene(spec([
        {"type": "reaction", "above": "HBr", "reagents": ["Br"],
         "from": ["CC=C1CCC1"], "to": ["CC1(Br)CCCC1"]}
    ]))
    assert scene.molecules == 2


def test_unbalanced_reaction_is_refused():
    with pytest.raises(SpecError) as exc:
        render_scene(spec([
            {"type": "reaction", "reagents": ["Br"],
             "from": ["CC=C1CCC1"], "to": ["CC1(C)CCCC1Br"]}
        ]))
    message = str(exc.value)
    assert "does not balance" in message
    assert "C: 6 on the left, 7 on the right" in message


def test_reaction_with_several_reactants():
    scene = render_scene(spec([
        {"type": "reaction", "above": "cat.",
         "from": ["CC(=O)O", "CCO"], "to": ["CC(=O)OCC", "O"]}
    ]))
    assert scene.molecules == 4


def test_reaction_needs_both_sides():
    with pytest.raises(SpecError) as exc:
        render_scene(spec([{"type": "reaction", "from": ["CCO"]}]))
    assert "from" in str(exc.value) and "to" in str(exc.value)


def test_overlap_warning():
    from chemdraft.scene import Box, Leaf, _overlap_warnings

    a = Leaf(kind="text", spec={})
    a.box = Box(0, 0, 100, 40)
    b = Leaf(kind="molecule", spec={})
    b.box = Box(80, 20, 100, 40)
    c = Leaf(kind="text", spec={})
    c.box = Box(300, 300, 100, 40)

    warnings = _overlap_warnings([a, b, c])
    assert len(warnings) == 1
    assert "overlap" in warnings[0]


