#!/usr/bin/env python3
"""Prepare plotted data from local result tables. No input data is uploaded."""
import argparse
import json
import os
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def local_path(value):
    if "://" in value:
        raise argparse.ArgumentTypeError("Use a localized file or directory, not a cloud URI")
    return Path(value)


def read_manifest(root=ROOT):
    value = json.loads((root / "preparation/manifest.json").read_text())
    return value.get("steps", []) if isinstance(value, dict) else value


def plan(steps, selectors):
    by_key = {x["key"]: x for x in steps}
    if len(by_key) != len(steps):
        raise ValueError("Duplicate preparation step keys")
    selected = []
    for selector in selectors:
        matches = [x["key"] for x in steps if x["key"] == selector or selector in x.get("figures", [])]
        if not matches:
            raise ValueError("Unknown preparation step or figure: " + selector)
        selected.extend(matches)
    result, done, visiting = [], set(), set()

    def visit(key):
        if key in visiting:
            raise ValueError("Cycle in preparation dependencies: " + key)
        if key in done:
            return
        if key not in by_key:
            raise ValueError("Unknown preparation dependency: " + key)
        visiting.add(key)
        for parent in by_key[key].get("dependencies", []):
            visit(parent)
        visiting.remove(key)
        done.add(key)
        result.append(by_key[key])

    for key in selected:
        visit(key)
    return result


def validate_manifest(steps, root=ROOT):
    plan(steps, [x["key"] for x in steps])
    for step in steps:
        script = Path(step["script"])
        if script.is_absolute() or ".." in script.parts or not (root / script).is_file():
            raise ValueError("Missing or unsafe preparation script: " + str(script))
        source = (root / script).read_text()
        if "/Users/" in source or "/private/tmp/" in source:
            raise ValueError("Personal path in preparation script: " + str(script))
        for output in step.get("outputs", []):
            path = Path(output["path"])
            if output["root"] not in {"paper", "source"} or path.is_absolute() or ".." in path.parts:
                raise ValueError("Unsafe preparation output: " + str(output))


def candidates(item, roots, env):
    root = item["root"]
    if item.get("environment_variable") and env.get(item["environment_variable"]):
        return [Path(env[item["environment_variable"]])]
    paths = [item["path"]] if item.get("path") else []
    paths += item.get("alternatives", [])
    if root.startswith("env:"):
        override = env.get(root[4:], "")
        if override:
            return [Path(override)]
        root = item.get("fallback_root", "input")
    result = []
    for relative in paths:
        if root in {"paper", "source"} and item.get("read_mode") != "original":
            result.append(roots["prepared"] / root / relative)
        path = roots[root] / relative
        result.extend(sorted(path.parent.glob(path.name)) if any(c in path.name for c in "*?[") else [path])
    return result


def check_inputs(steps, roots, env):
    produced = set()
    missing = []
    for step in steps:
        for item in step.get("inputs", []):
            if not item.get("required", True):
                continue
            outputs = [(item["root"], p) for p in [item.get("path", ""), *item.get("alternatives", [])]]
            if item.get("read_mode") != "original" and any(key in produced for key in outputs):
                continue
            if not any(p.is_file() and os.access(p, os.R_OK) for p in candidates(item, roots, env)):
                missing.append(step["key"] + ": " + item["root"] + "/" + (item.get("path") or ", ".join(item.get("alternatives", []))))
        produced.update((x["root"], x["path"]) for x in step.get("outputs", []))
    if missing:
        raise FileNotFoundError("Missing preparation inputs:\n" + "\n".join(missing))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("selectors", nargs="*", help="Figure names or preparation step keys")
    parser.add_argument("--list", action="store_true")
    parser.add_argument("--check-manifest", action="store_true")
    parser.add_argument("--check-inputs", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    for name in ["paper", "source", "input", "susie", "output"]:
        parser.add_argument("--" + name + "-root", type=local_path)
    parser.add_argument("--eas-eqtl-file", type=local_path)
    parser.add_argument("--topmed-cs-file", type=local_path)
    parser.add_argument("--protein-bed-file", type=local_path)
    parser.add_argument("--protein-covariates-file", type=local_path)
    parser.add_argument("--fig4-membership-file", type=local_path)
    parser.add_argument("--fig4-signal-file", type=local_path)
    args = parser.parse_args(argv)
    steps = read_manifest()
    validate_manifest(steps)
    if args.check_manifest:
        print(f"Validated {len(steps)} preparation steps.")
        return
    if args.list or not args.selectors:
        for step in steps:
            print(step["key"], ",".join(step["figures"]), step["description"], sep="\t")
        return
    selected = plan(steps, args.selectors)
    if args.dry_run:
        for step in selected:
            print(step["key"] + ": Rscript " + step["script"])
        return
    env = os.environ.copy()
    variables = {"paper": "QTL_PAPER_ROOT", "source": "QTL_SOURCE_ROOT", "input": "QTL_INPUT_ROOT",
                 "susie": "QTL_SUSIE_ROOT", "output": "QTL_PREP_ROOT"}
    for name, variable in variables.items():
        value = getattr(args, name + "_root")
        if value is not None:
            env[variable] = str(value.resolve())
    if args.eas_eqtl_file is not None:
        env["QTL_FIG1_EAS_EQTL_FILE"] = str(args.eas_eqtl_file.resolve())
    if args.topmed_cs_file is not None:
        env["QTL_TOPMED_CS_FILE"] = str(args.topmed_cs_file.resolve())
    for name, variable in {
        "protein_bed_file": "QTL_FIG1_PROTEIN_BED_FILE",
        "protein_covariates_file": "QTL_FIG1_PROTEIN_COVARIATES_FILE",
        "fig4_membership_file": "QTL_FIG4_MEMBERSHIP",
        "fig4_signal_file": "QTL_FIG4_SIGNAL",
    }.items():
        value = getattr(args, name)
        if value is not None:
            env[variable] = str(value.resolve())
    source = env.get("QTL_SOURCE_ROOT", str(ROOT))
    input_root = env.get("QTL_INPUT_ROOT", source)
    roots = {"source": Path(source), "input": Path(input_root),
             "paper": Path(env.get("QTL_PAPER_ROOT", str(ROOT))),
             "susie": Path(env.get("QTL_SUSIE_ROOT", input_root))}
    roots["prepared"] = Path(env.get("QTL_PREP_ROOT", str(ROOT / "notebook_outputs/prepared")))
    check_inputs(selected, roots, env)
    if args.check_inputs:
        print("All required preparation inputs are readable.")
        return
    for step in selected:
        print("Prepare " + step["key"], flush=True)
        expression = "source(commandArgs(trailingOnly=TRUE)[[1]]); prep_finish(commandArgs(trailingOnly=TRUE)[[2]])"
        subprocess.run(["Rscript", "-e", expression, step["script"], step["key"]], cwd=ROOT, env=env, check=True)


if __name__ == "__main__":
    main()
