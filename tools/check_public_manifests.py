#!/usr/bin/env python3
"""Reject concrete GCS addresses in files prepared for publication."""

from pathlib import Path
import re
import subprocess
import sys


GCS_ADDRESS = re.compile(r"gs://([A-Za-z0-9][A-Za-z0-9._-]*)(?:/[^\s\"'`<>)]*)?")


def check_file(path: Path, root: Path) -> list[str]:
    """Return file/line diagnostics without disclosing matched addresses."""
    try:
        content = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        return []
    relative = path.relative_to(root)
    synthetic_test = relative.parts[0] == "tests" and path.suffix in {".R", ".py"}
    issues = []
    for number, line in enumerate(content.splitlines(), 1):
        matches = GCS_ADDRESS.finditer(line)
        if any(not (synthetic_test and match.group(1) == "bucket") for match in matches):
            issues.append(f"{relative.as_posix()}:{number}: remove the GCS address from this public file")
    return issues


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    names = subprocess.check_output(
        ["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z"], cwd=root,
    ).decode("utf-8").split("\0")
    issues = []
    checked = 0
    for name in sorted(set(filter(None, names))):
        path = root / name
        if path.is_file():
            checked += 1
            issues.extend(check_file(path, root))
    if issues:
        print("\n".join(issues), file=sys.stderr)
        print(f"Public-file check failed: {len(issues)} line(s) contain a GCS address.", file=sys.stderr)
        return 1
    print(f"Public-file check passed: {checked} files checked; no GCS addresses.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
