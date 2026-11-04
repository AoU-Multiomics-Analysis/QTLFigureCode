#!/usr/bin/env python3
"""Check current notebook code and outputs without access to controlled inputs."""
import argparse
import base64
import hashlib
import json
from pathlib import Path
import re
import shutil
import tempfile

from build_current_notebooks import ROOT, build, read_index


def check(require_executed=False):
    failures = []
    for key, entry in read_index().items():
        path = ROOT / entry["path"]
        try:
            nb = json.loads(path.read_text())
            assert nb["nbformat"] == 4
            assert nb["metadata"]["qtl_reference_date"] == "2026-11-02"
            assert nb["metadata"]["kernelspec"]["name"] == "ir"
            ids = [c["id"] for c in nb["cells"]]
            assert len(ids) == len(set(ids)), "Duplicate cell IDs"
            assert all(isinstance(c["metadata"], dict) for c in nb["cells"])
            with tempfile.TemporaryDirectory() as directory:
                temp = Path(directory)
                for name in ["notebooks/current/00_setup.R", entry["script"]]:
                    target = temp / name
                    target.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(ROOT / name, target)
                expected_path, _ = build(key, entry, root=temp)
                expected = json.loads(expected_path.read_text())
                actual_code = [c["source"] for c in nb["cells"] if c["cell_type"] == "code"]
                expected_code = [c["source"] for c in expected["cells"] if c["cell_type"] == "code"]
                assert actual_code == expected_code, "Notebook code differs from its source script"
            images = []
            for c in nb["cells"]:
                source = "".join(c["source"])
                assert "/Users/" not in source and "/private/tmp/" not in source, "Personal path in cell"
                if c["cell_type"] != "code":
                    continue
                if require_executed:
                    assert isinstance(c["execution_count"], int), "Unexecuted code cell"
                for output in c["outputs"]:
                    assert output["output_type"] in {"stream", "display_data"}, "Unexpected or failed cell output"
                    if output["output_type"] == "display_data":
                        assert set(output["data"]) == {"image/png"}, "Only aggregate plot images are permitted"
                        assert isinstance(output["metadata"], dict)
                        image = base64.b64decode(output["data"]["image/png"])
                        assert image.startswith(b"\x89PNG\r\n\x1a\n"), "Invalid PNG output"
                        images.append(image)
                    else:
                        text = "".join(output["text"]) if isinstance(output["text"], list) else output["text"]
                        assert not re.search(r"\b(research_id|participant_id|person_id)\b", text), "Participant fields in output"
                        assert "/Users/" not in text and "/private/tmp/" not in text, "Personal path in output"
                        assert "conversion failure" not in text and "invalid font" not in text, "Plot text failed to render"
            if require_executed:
                execution = nb["metadata"]["qtl_execution"]
                assert execution["source_inputs_unchanged"] is True
                assert images, "No rendered plot"
                if entry.get("reference_sha256"):
                    assert hashlib.sha256(images[-1]).hexdigest() == entry["reference_sha256"], "Submission reference image differs"
            print(f"PASS {key}: {len(actual_code)} code cells, {len(images)} images")
        except (AssertionError, KeyError, ValueError, OSError) as error:
            failures.append(f"{key}: {error}")
    if failures:
        raise SystemExit("\n".join(failures))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--require-executed", action="store_true")
    check(parser.parse_args().require_executed)
