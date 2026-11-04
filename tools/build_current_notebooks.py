#!/usr/bin/env python3
"""Build current R notebooks and keep outputs for code that has not changed."""
import argparse
import hashlib
import json
import os
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read_index(root=ROOT):
    return json.loads((root / "notebooks/current/index.json").read_text())


def expand_source(name, root=ROOT, stack=()):
    """Expand explicit repository-local R includes before splitting cells."""
    root = Path(root).resolve()
    path = (root / name).resolve()
    if not path.is_relative_to(root):
        raise ValueError("Notebook include is outside the repository: " + str(name))
    relative = path.relative_to(root).as_posix()
    if relative in stack:
        raise ValueError("Notebook include cycle: " + " -> ".join((*stack, relative)))
    lines, dependencies = [], [relative]
    for line in path.read_text().splitlines():
        match = re.fullmatch(r"# @include (.+?)\s*", line)
        if match:
            included, files = expand_source(match[1], root, (*stack, relative))
            lines.append(included)
            dependencies.extend(files)
        else:
            lines.append(line)
    return "\n".join(lines) + "\n", list(dict.fromkeys(dependencies))


def source_files(entry, root=ROOT):
    files = []
    for name in ["notebooks/current/00_setup.R", entry["script"]]:
        _, dependencies = expand_source(name, root)
        files.extend(dependencies)
    return list(dict.fromkeys(files))


def cell(kind, text, label):
    text = text.strip() + "\n"
    result = {"cell_type": kind, "id": hashlib.sha256((kind + label).encode()).hexdigest()[:16],
              "metadata": {}, "source": text.splitlines(keepends=True)}
    if kind == "code":
        result.update(execution_count=None, outputs=[])
    return result


def script_sections(text):
    """Read authored R sections. Lines beginning '#|' become method notes."""
    title, notes, lines = "Plot code", [], []
    for line in text.splitlines():
        if line.startswith("# %% "):
            if any(s.strip() for s in lines):
                yield title, "\n".join(notes), "\n".join(lines)
            title, notes, lines = line[5:].strip(), [], []
        elif line.startswith("#| "):
            notes.append(line[3:])
        else:
            lines.append(line)
    if any(s.strip() for s in lines):
        yield title, "\n".join(notes), "\n".join(lines)


def build(key, entry, root=ROOT, clear_outputs=False):
    root = Path(root)
    script, script_files = expand_source(entry["script"], root)
    setup, setup_files = expand_source("notebooks/current/00_setup.R", root)
    cells = [cell("markdown", f"# {entry['title']}\n\n{entry['overview']}", "intro")]
    cells.append(cell("markdown", """## Run the notebook

Use the R kernel. Start in the repository root. Set `QTL_INPUT_ROOT` to the
source result tables and `QTL_SUSIE_ROOT` to the separate fine-mapping tables.
Set `QTL_PAPER_ROOT` to the local reference directory and `QTL_SOURCE_ROOT`
to the figure checkout. Set these variables in R before running the cells.
Paths can contain spaces.

Run the cells from top to bottom. The cells show input reads, filters, joins,
counts, models, and plot construction. Intermediate objects stay in this R
session. New PDFs and PNGs go to
`notebook_outputs/`. The plotting scripts are in `notebooks/current/`.
""", "run"))
    if entry.get("preparation"):
        notebook_dir = (root / entry["path"]).parent
        links = []
        for step in entry["preparation"]:
            target = os.path.relpath(root / step["script"], notebook_dir)
            links.append(f"- [{step['key']}]({target}): {step['description']}")
        instructions = (
            "## Source scripts\n\n"
            "The processing steps below follow these preparation scripts. "
            "The source data are not included in Git.\n\n"
            + "\n".join(links) + "\n\n"
            "The notebook runs the processing in its own cells. "
            "The preparation guide lists the original input requirements and missing inputs."
        )
        cells.append(cell("markdown", instructions, "preparation"))
    cells += [cell("markdown", "## Setup\n\nValidate the local input paths and load the plot libraries.", "setup-note"),
              cell("code", setup, "setup")]
    for index, (title, notes, body) in enumerate(script_sections(script)):
        label = f"section-{index}-{title}"
        cells += [cell("markdown", f"## {title}\n\n{notes}", label + "-note"),
                  cell("code", body, label)]
    out = root / entry["path"]
    execution = None
    if out.exists() and not clear_outputs:
        old = json.loads(out.read_text())
        old_code = [c for c in old["cells"] if c["cell_type"] == "code"]
        position = 0
        unchanged = True
        for c in cells:
            if c["cell_type"] != "code":
                continue
            previous = old_code[position] if position < len(old_code) else None
            unchanged = unchanged and previous is not None and c["source"] == previous["source"]
            if unchanged:
                c["outputs"] = previous.get("outputs", [])
                c["execution_count"] = previous.get("execution_count")
            position += 1
        if unchanged and position == len(old_code):
            execution = old.get("metadata", {}).get("qtl_execution")
    nb = {"cells": cells, "metadata": {
        "kernelspec": {"name": "ir", "display_name": "R", "language": "R"},
        "language_info": {"name": "R", "file_extension": ".r", "mimetype": "text/x-r-source"},
        "qtl_reference_date": "2026-11-02", "qtl_figure": key,
        "qtl_source_script": entry["script"],
        "qtl_source_files": list(dict.fromkeys(setup_files + script_files))}, "nbformat": 4, "nbformat_minor": 5}
    if execution is not None:
        nb["metadata"]["qtl_execution"] = execution
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(nb, indent=1) + "\n")
    return out, len(cells)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("figures", nargs="*")
    parser.add_argument("--clear-outputs", action="store_true")
    parser.add_argument("--list", action="store_true")
    args = parser.parse_args(argv)
    index = read_index()
    if args.list:
        for key, entry in index.items():
            print(key, entry["path"])
        return
    targets = args.figures or list(index)
    unknown = set(targets) - set(index)
    if unknown:
        parser.error("Unknown figure(s): " + ", ".join(sorted(unknown)))
    for key in targets:
        path, count = build(key, index[key], clear_outputs=args.clear_outputs)
        print(f"Wrote {path.relative_to(ROOT)} ({count} cells)")


if __name__ == "__main__":
    main()
