"""Shared configuration for the AlphaGenome scoring/prediction scripts.

All external paths and the API key live here; override any with an environment
variable of the same name. Run with the interpreter that has `alphagenome`
installed (see README): e.g. /Users/evinmpadhi/borzoi_run/venv/bin/python.
"""
import os

# AlphaGenome results directory (parquets, tsvs, figures)
AG_DIR = os.environ.get(
    "FIG5_AG", "/Users/evinmpadhi/Fig5_QTL_paper/alphagenome")
os.makedirs(AG_DIR, exist_ok=True)

# AlphaGenome API key — set AG_KEY in the environment; do NOT hard-code secrets.
AG_KEY = os.environ.get("AG_KEY")


def require_key() -> str:
    """Return the API key or raise a clear error if it is not set."""
    if not AG_KEY:
        raise SystemExit(
            "AG_KEY environment variable is not set. Export your AlphaGenome "
            "API key, e.g.  export AG_KEY=...  before running this script.")
    return AG_KEY


def path(*parts: str) -> str:
    """Path helper rooted at AG_DIR."""
    return os.path.join(AG_DIR, *parts)
