"""chemdraft command line.

One job: turn a SMILES string into a figure that provably shows that molecule.
Writes nothing when it cannot prove it.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .chemistry import ChemError
from .render import DEFAULT_BOND_PX, RenderRequest, render


def _emit(payload: dict, code: int) -> int:
    json.dump(payload, sys.stdout, indent=2)
    sys.stdout.write("\n")
    return code


def cmd_render(args: argparse.Namespace) -> int:
    from .scene import SpecError, load_spec, render_scene

    try:
        spec = load_spec(args.spec)
        if args.no_reference:
            spec.setdefault("style", {})["reference"] = False
        scene = render_scene(spec)
    except SpecError as exc:
        return _emit({"ok": False, "error": str(exc), "written": []}, 2)
    except Exception as exc:  # noqa: BLE001 - fail closed
        return _emit(
            {"ok": False,
             "error": f"internal error, nothing written: {type(exc).__name__}: {exc}",
             "written": []},
            4,
        )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_bytes(scene.png)
    written = [str(out)]
    svg_path = None
    if args.svg:
        svg_path = Path(args.svg) if isinstance(args.svg, str) else out.with_suffix(".svg")
        svg_path.write_text(scene.svg, encoding="utf-8")
        written.append(str(svg_path))

    return _emit(
        {
            "ok": True,
            "molecules": scene.molecules,
            "png": str(out),
            "svg": str(svg_path) if svg_path else None,
            "width": scene.width,
            "height": scene.height,
            "caveats": scene.caveats,
            "written": written,
        },
        0,
    )


def cmd_smiles(args: argparse.Namespace) -> int:
    out = Path(args.out) if args.out else Path("structure.png")
    try:
        result = render(
            RenderRequest(
                smiles=args.smiles,
                mode=args.mode,
                highlight=args.highlight,
                bond_px=args.bond_px,
                color=args.color,
                reference=not args.no_reference,
            )
        )
    except ChemError as exc:
        return _emit({"ok": False, "error": str(exc), "written": []}, 2)
    except Exception as exc:  # noqa: BLE001 - fail closed, never write
        return _emit(
            {
                "ok": False,
                "error": f"internal error, nothing written: {type(exc).__name__}: {exc}",
                "written": [],
            },
            4,
        )

    if not result.report.ok:
        # Do not write a figure we cannot vouch for.
        return _emit(
            {
                "ok": False,
                "error": "the drawing did not pass verification, no file written",
                "problems": result.report.problems,
                "checks_passed": result.report.checks_run,
                "written": [],
            },
            3,
        )

    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_bytes(result.png)
    written = [str(out)]
    svg_path = None
    if args.svg:
        svg_path = Path(args.svg) if isinstance(args.svg, str) else out.with_suffix(".svg")
        svg_path.write_text(result.svg, encoding="utf-8")
        written.append(str(svg_path))

    d = result.description
    return _emit(
        {
            "ok": True,
            "smiles": d.smiles,
            "canonical": d.canonical,
            "formula": d.formula,
            "mw": d.mw,
            "charge": d.charge,
            "atoms": d.atoms,
            "bonds": d.bonds,
            "stereocentres": d.stereocentres,
            "unassigned_stereocentres": d.unassigned_stereocentres,
            "double_bond_stereo": d.double_bond_stereo,
            "png": str(out),
            "svg": str(svg_path) if svg_path else None,
            "width": result.width,
            "height": result.height,
            "depiction": result.depiction,
            "embed": f"![[{out.name}]]",
            "verified": True,
            "checks": result.report.checks_run,
            "warnings": result.report.problems + result.report.notes + result.notes,
            "written": written,
        },
        0,
    )


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="chemdraft",
        description="Draw a molecule from SMILES, verified before it is written.",
    )
    sub = p.add_subparsers(dest="command", required=True)

    s = sub.add_parser("smiles", help="draw one molecule")
    s.add_argument("smiles", help="SMILES string, quoted")
    s.add_argument("-o", "--out", help="PNG path (default structure.png)")
    s.add_argument("--svg", nargs="?", const=True, help="also write SVG")
    s.add_argument("--mode", default="skeletal", choices=["skeletal", "lewis"],
                   help="skeletal hides carbons, lewis shows every atom")
    s.add_argument("--highlight", help="SMARTS of the part to highlight")
    s.add_argument("--bond-px", type=int, default=DEFAULT_BOND_PX)
    s.add_argument("--no-reference", action="store_true",
                   help="do not look up a published depiction, lay out automatically")
    s.add_argument("--color", action="store_true",
                   help="colour heteroatoms instead of black and white")
    s.set_defaults(func=cmd_smiles)

    r = sub.add_parser("render", help="draw a scene from a JSON spec")
    r.add_argument("spec", help="path to a JSON scene description")
    r.add_argument("-o", "--out", default="scene.png", help="PNG path")
    r.add_argument("--svg", nargs="?", const=True, help="also write SVG")
    r.add_argument("--no-reference", action="store_true",
                   help="do not look up published depictions")
    r.set_defaults(func=cmd_render)

    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
