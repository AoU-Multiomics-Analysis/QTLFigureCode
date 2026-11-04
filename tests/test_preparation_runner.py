"""Check preparation planning without controlled source data."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("preparation", ROOT / "tools/prepare_figure_data.py")
preparation = importlib.util.module_from_spec(spec)
spec.loader.exec_module(preparation)


class PreparationRunnerTest(unittest.TestCase):
    def test_order_dependencies_once(self):
        steps = [dict(key="filter", figures=["Fig1"], dependencies=[]),
                 dict(key="count", figures=["Fig1"], dependencies=["filter"]),
                 dict(key="plot_input", figures=["S11"], dependencies=["count", "filter"])]
        self.assertEqual([x["key"] for x in preparation.plan(steps, ["plot_input"])],
                         ["filter", "count", "plot_input"])

    def test_select_figure_and_reject_unknown_steps(self):
        steps = [dict(key="filter", figures=["Fig1"], dependencies=[])]
        self.assertEqual(len(preparation.plan(steps, ["Fig1"])), 1)
        with self.assertRaisesRegex(ValueError, "Unknown"):
            preparation.plan(steps, ["absent"])

    def test_cycle_stops_before_execution(self):
        steps = [dict(key="a", figures=["Fig1"], dependencies=["b"]),
                 dict(key="b", figures=["Fig1"], dependencies=["a"])]
        with self.assertRaisesRegex(ValueError, "[Cc]ycle"):
            preparation.plan(steps, ["a"])

    def test_missing_script_or_personal_path_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = [dict(key="a", script="preparation/missing.R", figures=["Fig1"], dependencies=[], inputs=[], outputs=[])]
            with self.assertRaisesRegex(ValueError, "script"):
                preparation.validate_manifest(manifest, root)
            (root / "preparation").mkdir()
            (root / "preparation/missing.R").write_text('root <- "/Users/a/private"\n')
            with self.assertRaisesRegex(ValueError, "[Pp]ersonal"):
                preparation.validate_manifest(manifest, root)

    def test_missing_input_is_reported_before_subprocess(self):
        with tempfile.TemporaryDirectory(prefix="qtl source roots ") as directory:
            roots = {"input": Path(directory), "paper": Path(directory),
                     "source": Path(directory), "susie": Path(directory),
                     "prepared": Path(directory) / "prepared"}
            step = dict(key="filter", inputs=[dict(root="input", path="missing.tsv")], outputs=[])
            with self.assertRaisesRegex(FileNotFoundError, "missing.tsv"):
                preparation.check_inputs([step], roots, {})
            (Path(directory) / "missing.tsv").write_text("fixture\n")
            preparation.check_inputs([step], roots, {})

    def test_original_reference_cannot_be_satisfied_by_overlay_or_planned_output(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            roots = {name: root / name for name in ["paper", "source", "input", "susie", "prepared"]}
            generated = roots["prepared"] / "paper/template.rds"
            generated.parent.mkdir(parents=True)
            generated.write_text("prior output")
            reference = dict(key="connect", inputs=[dict(root="paper", path="template.rds", read_mode="original")], outputs=[])
            with self.assertRaisesRegex(FileNotFoundError, "template.rds"):
                preparation.check_inputs([reference], roots, {})
            producer = dict(key="producer", inputs=[], outputs=[dict(root="paper", path="template.rds")])
            with self.assertRaisesRegex(FileNotFoundError, "template.rds"):
                preparation.check_inputs([producer, reference], roots, {})
            roots["paper"].mkdir()
            (roots["paper"] / "template.rds").write_text("original")
            preparation.check_inputs([producer, reference], roots, {})


if __name__ == "__main__":
    unittest.main()
